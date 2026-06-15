#!/usr/bin/env node

import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';

const BASES = ['A', 'C', 'G', 'T', 'N'];

function main() {
  const configPath = path.resolve(process.argv[2] || 'benchmarks/public/crispresso2_base_editor/config.json');
  const configDir = path.dirname(configPath);
  const config = JSON.parse(fs.readFileSync(configPath, 'utf8'));
  const resultDir = path.resolve(configDir, config.outputs?.resultDir || 'results');
  fs.mkdirSync(resultDir, { recursive: true });

  const reference = normalizeSeq(config.reference);
  if (reference.length < 20) throw new Error('Reference amplicon sequence must be at least 20 bp.');

  const settings = normalizeSettings(config.settings);
  const assay = deriveAssayTargets(reference, config.assay || {}, config.targetPositions || '');
  const inputFile = path.resolve(configDir, config.dataset?.file || config.fastq || '');
  if (!fs.existsSync(inputFile)) {
    throw new Error(`FASTQ file not found: ${inputFile}`);
  }

  const sampleName = config.dataset?.sampleName || sampleNameFromFile(inputFile);
  const inputFileLabel = config.dataset?.file || path.relative(configDir, inputFile).replace(/\\/g, '/');
  const fastqText = readFastq(inputFile);
  const parsed = parseFastqText(fastqText, settings);
  const sample = analyzeSample({
    sample: sampleName,
    files: [path.basename(inputFile)],
    roles: ['R1']
  }, parsed, reference, settings);

  const run = buildRunResult([sample], reference, settings, assay, config, inputFileLabel);
  writeOutputs(resultDir, run);
  printSummary(run, resultDir);
}

function normalizeSettings(settings = {}) {
  return {
    minQ: Math.max(0, Number(settings.minQ ?? 20) || 0),
    minLen: Math.max(1, Number(settings.minLen ?? 30) || 1),
    readLimit: Math.max(1, Number(settings.readLimit ?? 50000) || 1),
    minIdentity: normalizeMinIdentity(settings.minIdentity ?? 0.8),
    dropN: settings.dropN !== false,
    minDepth: Math.max(0, Number(settings.minDepth ?? 0) || 0),
    signalMode: ['all', 'substitution', 'indel'].includes(settings.signalMode) ? settings.signalMode : 'all'
  };
}

function normalizeMinIdentity(value) {
  const n = Number(value);
  if (!Number.isFinite(n)) return 0.8;
  return Math.max(0, Math.min(1, n > 1 ? n / 100 : n));
}

function readFastq(filePath) {
  const data = fs.readFileSync(filePath);
  if (filePath.toLowerCase().endsWith('.gz')) return zlib.gunzipSync(data).toString('utf8');
  return data.toString('utf8');
}

function sampleNameFromFile(filePath) {
  return path.basename(filePath).replace(/(\.fastq|\.fq)(\.gz)?$/i, '');
}

function normalizeSeq(value) {
  return String(value || '').toUpperCase().replace(/[^ACGTN]/g, '');
}

function normalizePattern(value) {
  return String(value || '').toUpperCase().replace(/[^ACGTRYSWKMBDHVN]/g, 'N');
}

function parseTargets(value) {
  const targets = new Set();
  String(value || '').split(/[,\s;]+/).filter(Boolean).forEach((part) => {
    const range = part.match(/^(\d+)-(\d+)$/);
    if (range) {
      const start = Number(range[1]);
      const end = Number(range[2]);
      for (let pos = Math.min(start, end); pos <= Math.max(start, end); pos += 1) targets.add(pos);
    } else {
      const pos = Number(part.replace(/\D/g, ''));
      if (pos) targets.add(pos);
    }
  });
  return targets;
}

function deriveAssayTargets(reference, assay, manualTargets) {
  const manual = parseTargets(manualTargets);
  const fallback = (warning = '') => ({
    targets: manual,
    source: manual.size ? 'manual / fallback target positions' : 'no target positions',
    warning,
    details: {}
  });

  if (assay.type === 'cas') {
    const spacer = normalizeSeq(assay.spacer || '');
    const pamPattern = normalizePattern(assay.pam || 'NGG');
    const spacerWindow = parseTargets(assay.spacerWindow || '4-8');
    if (!spacer) return fallback('Cas mode: no spacer/guide sequence was provided.');
    const match = findSpacerInReference(reference, spacer, pamPattern);
    if (!match) return fallback('Cas mode: spacer was not found in the reference.');
    const relativePositions = spacerWindow.size ? spacerWindow : parseTargets(`1-${spacer.length}`);
    const targets = new Set();
    relativePositions.forEach((pos) => {
      if (pos < 1 || pos > spacer.length) return;
      const refPosition = match.orientation === 'forward'
        ? match.start + pos - 1
        : match.start + spacer.length - pos;
      if (refPosition >= 1 && refPosition <= reference.length) targets.add(refPosition);
    });
    return {
      targets,
      source: `Cas ${match.orientation} spacer match - PAM ${match.pamMatched ? 'matched' : 'not confirmed'} - spacer window ${assay.spacerWindow || '4-8'}`,
      warning: match.pamMatched ? '' : 'Cas mode: spacer was found, but the PAM pattern was not confirmed.',
      details: {
        spacer,
        pamPattern,
        spacerWindow: assay.spacerWindow || '4-8',
        orientation: match.orientation,
        spacerStart: match.start,
        spacerEnd: match.end,
        expectedEdit: assay.expectedEdit || null
      }
    };
  }

  return fallback('');
}

function findSpacerInReference(reference, spacer, pamPattern) {
  const candidates = [];
  collectSpacerCandidates(candidates, reference, spacer, pamPattern, 'forward');
  collectSpacerCandidates(candidates, reference, reverseComplement(spacer), pamPattern, 'reverse');
  candidates.sort((a, b) => b.score - a.score || Math.abs(centerOf(reference) - centerOfMatch(a)) - Math.abs(centerOf(reference) - centerOfMatch(b)));
  return candidates[0] || null;
}

function collectSpacerCandidates(candidates, reference, query, pamPattern, orientation) {
  if (!query) return;
  let start = reference.indexOf(query);
  while (start !== -1) {
    const oneBasedStart = start + 1;
    const pamInfo = scorePam(reference, start, query.length, pamPattern, orientation);
    candidates.push({
      start: oneBasedStart,
      end: oneBasedStart + query.length - 1,
      orientation,
      pamMatched: pamInfo.matched,
      score: query.length + pamInfo.score
    });
    start = reference.indexOf(query, start + 1);
  }
}

function scorePam(reference, startZero, spacerLength, pamPattern, orientation) {
  if (!pamPattern) return { matched: true, score: 0 };
  const pamLength = pamPattern.length;
  let pamSeq = '';
  let pattern = pamPattern;
  if (orientation === 'forward') {
    pamSeq = reference.slice(startZero + spacerLength, startZero + spacerLength + pamLength);
  } else {
    pamSeq = reference.slice(Math.max(0, startZero - pamLength), startZero);
    pattern = reverseComplementPattern(pamPattern);
  }
  const matched = pamSeq.length === pamLength && matchIupac(pamSeq, pattern);
  return { matched, score: matched ? 100 : -8 };
}

function centerOf(reference) {
  return reference.length / 2;
}

function centerOfMatch(match) {
  return (match.start + match.end) / 2;
}

function matchIupac(sequence, pattern) {
  const map = {
    A: 'A', C: 'C', G: 'G', T: 'T', N: 'ACGT',
    R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC',
    B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG'
  };
  if (sequence.length !== pattern.length) return false;
  for (let i = 0; i < pattern.length; i += 1) {
    const allowed = map[pattern[i]] || pattern[i];
    if (!allowed.includes(sequence[i])) return false;
  }
  return true;
}

function reverseComplementPattern(pattern) {
  const map = { A: 'T', C: 'G', G: 'C', T: 'A', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D', N: 'N' };
  return String(pattern || '').split('').reverse().map((base) => map[base] || 'N').join('');
}

function parseFastqText(text, settings) {
  const lines = String(text || '').replace(/\r\n/g, '\n').replace(/\r/g, '\n').split('\n');
  const estimatedRecords = Math.floor(lines.length / 4);
  const limit = Math.min(settings.readLimit, estimatedRecords || settings.readLimit);
  const unique = new Map();
  const readLengths = [];
  let totalRecords = 0;
  let keptRecords = 0;
  let malformedRecords = 0;
  let droppedShort = 0;
  let droppedQuality = 0;
  let droppedN = 0;
  let totalQuality = 0;
  let qualityReads = 0;
  let basesKept = 0;

  for (let i = 0; i + 3 < lines.length && totalRecords < limit; i += 4) {
    const header = lines[i];
    const seqLine = lines[i + 1];
    const plus = lines[i + 2];
    const qualLine = lines[i + 3];
    if (!header || header[0] !== '@' || !plus || plus[0] !== '+') {
      malformedRecords += 1;
      continue;
    }
    totalRecords += 1;
    const seq = normalizeSeq(seqLine);
    const qual = String(qualLine || '').trim();
    if (!seq || qual.length < seq.length) {
      malformedRecords += 1;
      continue;
    }
    if (seq.length < settings.minLen) {
      droppedShort += 1;
      continue;
    }
    if (settings.dropN && seq.includes('N')) {
      droppedN += 1;
      continue;
    }
    const meanQ = meanQuality(qual.slice(0, seq.length));
    totalQuality += meanQ;
    qualityReads += 1;
    if (meanQ < settings.minQ) {
      droppedQuality += 1;
      continue;
    }
    keptRecords += 1;
    basesKept += seq.length;
    readLengths.push(seq.length);
    const current = unique.get(seq) || { sequence: seq, count: 0, qualitySum: 0 };
    current.count += 1;
    current.qualitySum += meanQ;
    unique.set(seq, current);
  }

  return {
    estimatedRecords,
    readLimitReached: estimatedRecords > limit,
    totalRecords,
    keptRecords,
    malformedRecords,
    droppedShort,
    droppedQuality,
    droppedN,
    unique,
    uniqueReads: unique.size,
    meanQuality: qualityReads ? totalQuality / qualityReads : 0,
    meanReadLength: keptRecords ? basesKept / keptRecords : 0,
    medianReadLength: median(readLengths)
  };
}

function meanQuality(qual) {
  let sum = 0;
  for (let i = 0; i < qual.length; i += 1) sum += Math.max(0, qual.charCodeAt(i) - 33);
  return qual.length ? sum / qual.length : 0;
}

function median(values) {
  if (!values.length) return 0;
  const copy = values.slice().sort((a, b) => a - b);
  return copy[Math.floor(copy.length / 2)];
}

function analyzeSample(group, parsedGroup, reference, settings) {
  const stats = Array.from({ length: reference.length }, (_, idx) => ({
    position: idx + 1,
    refBase: reference[idx],
    substitutionReads: 0,
    deletionReads: 0,
    insertionReads: 0,
    editedReads: 0,
    coveredReads: 0,
    A: 0,
    C: 0,
    G: 0,
    T: 0,
    N: 0
  }));
  const alleles = [];
  let alignedReads = 0;
  let noIndelReads = 0;
  let editedReads = 0;
  let lowIdentityReads = 0;
  let totalIdentity = 0;
  let identityWeight = 0;
  const uniqueRows = Array.from(parsedGroup.unique.values()).sort((a, b) => b.count - a.count);

  for (const row of uniqueRows) {
    const alignment = alignRead(reference, row.sequence);
    if (alignment.identity < settings.minIdentity) {
      lowIdentityReads += row.count;
      continue;
    }
    const edits = collectEdits(alignment);
    alignedReads += row.count;
    totalIdentity += alignment.identity * row.count;
    identityWeight += row.count;
    if (!edits.hasIndel) noIndelReads += row.count;
    if (edits.hasAnyEdit) editedReads += row.count;

    if (edits.coveredStart && edits.coveredEnd) {
      for (let position = edits.coveredStart; position <= edits.coveredEnd; position += 1) {
        const stat = stats[position - 1];
        if (stat) stat.coveredReads += row.count;
      }
    }

    const editedPositions = new Set();
    for (const edit of edits.items) {
      const stat = stats[edit.position - 1];
      if (!stat) continue;
      if (edit.type === 'substitution') {
        stat.substitutionReads += row.count;
        const base = BASES.includes(edit.alt) ? edit.alt : 'N';
        stat[base] += row.count;
      } else if (edit.type === 'deletion') {
        stat.deletionReads += row.count;
      } else if (edit.type === 'insertion') {
        stat.insertionReads += row.count;
      }
      editedPositions.add(edit.position);
    }
    for (const position of editedPositions) {
      const stat = stats[position - 1];
      if (stat) stat.editedReads += row.count;
    }

    alleles.push({
      sample: group.sample,
      sequence: row.sequence,
      nReads: row.count,
      ratio: 0,
      orientation: alignment.orientation,
      identity: alignment.identity,
      editSignature: edits.signature || 'WT',
      mismatchCount: edits.mismatchCount,
      insertionCount: edits.insertionCount,
      deletionCount: edits.deletionCount,
      hasIndel: edits.hasIndel
    });
  }
  for (const allele of alleles) allele.ratio = alignedReads ? allele.nReads / alignedReads : 0;

  const rows = stats.map((stat) => {
    const displayReads = getSignalReads(stat, settings.signalMode);
    return {
      sample: group.sample,
      label: group.sample,
      position: stat.position,
      refBase: stat.refBase,
      substitutionReads: stat.substitutionReads,
      deletionReads: stat.deletionReads,
      insertionReads: stat.insertionReads,
      editedReads: stat.editedReads,
      displayReads,
      denominator: stat.coveredReads,
      percent: stat.coveredReads ? displayReads / stat.coveredReads * 100 : 0,
      coveredReads: stat.coveredReads,
      A: stat.A,
      C: stat.C,
      G: stat.G,
      T: stat.T,
      N: stat.N
    };
  });

  const maxPercent = Math.max(0, ...rows.map((row) => row.percent));
  const warnings = [];
  if (parsedGroup.keptRecords < settings.minDepth) warnings.push('below minimum reads');
  if (parsedGroup.malformedRecords) warnings.push('malformed FASTQ records');
  if (parsedGroup.readLimitReached) warnings.push('read limit applied');
  if (!parsedGroup.keptRecords) warnings.push('no QC-passed reads');
  if (lowIdentityReads) warnings.push('low-identity reads excluded');
  if (parsedGroup.keptRecords && !alignedReads) warnings.push('no reads passed alignment identity');

  return {
    sample: group.sample,
    files: group.files,
    roles: group.roles,
    qc: {
      estimatedRecords: parsedGroup.estimatedRecords,
      parsedRecords: parsedGroup.totalRecords,
      passedReads: parsedGroup.keptRecords,
      alignedReads,
      lowIdentityReads,
      noIndelReads,
      editedReads,
      uniqueReads: alleles.length,
      meanQuality: parsedGroup.meanQuality,
      medianReadLength: parsedGroup.medianReadLength,
      droppedShort: parsedGroup.droppedShort,
      droppedQuality: parsedGroup.droppedQuality,
      droppedN: parsedGroup.droppedN,
      malformedRecords: parsedGroup.malformedRecords,
      readLimitReached: parsedGroup.readLimitReached,
      meanIdentity: identityWeight ? totalIdentity / identityWeight : 0
    },
    positionRows: rows,
    alleles: alleles.sort((a, b) => b.nReads - a.nReads),
    maxPercent,
    warnings
  };
}

function getSignalReads(stat, mode) {
  if (mode === 'substitution') return stat.substitutionReads;
  if (mode === 'indel') return stat.deletionReads + stat.insertionReads;
  return stat.editedReads;
}

function alignRead(reference, sequence) {
  const forward = alignCandidate(reference, sequence, 'forward');
  const rcSeq = reverseComplement(sequence);
  const reverse = alignCandidate(reference, rcSeq, 'reverse');
  return reverse.score > forward.score ? reverse : forward;
}

function alignCandidate(reference, sequence, orientation) {
  const direct = directAlignment(reference, sequence, orientation);
  if (direct && direct.identity >= 0.86) return direct;
  return needlemanSemiGlobal(reference, trimLongRead(reference, sequence), orientation);
}

function trimLongRead(reference, sequence) {
  if (sequence.length <= reference.length + 80) return sequence;
  const window = Math.min(sequence.length, reference.length + 40);
  let bestStart = 0;
  let bestScore = -Infinity;
  const step = Math.max(1, Math.floor((sequence.length - window) / 80));
  for (let start = 0; start + window <= sequence.length; start += step) {
    const fragment = sequence.slice(start, start + window);
    let matches = 0;
    const limit = Math.min(reference.length, fragment.length);
    for (let i = 0; i < limit; i += 1) if (reference[i] === fragment[i]) matches += 1;
    if (matches > bestScore) {
      bestScore = matches;
      bestStart = start;
    }
  }
  return sequence.slice(bestStart, bestStart + window);
}

function directAlignment(reference, sequence, orientation) {
  if (reference.length !== sequence.length) return null;
  let matches = 0;
  let score = 0;
  for (let i = 0; i < reference.length; i += 1) {
    if (reference[i] === sequence[i]) {
      matches += 1;
      score += 2;
    } else {
      score -= 2;
    }
  }
  return {
    alignedRef: reference,
    alignedRead: sequence,
    score,
    orientation,
    refStartOffset: 0,
    identity: reference.length ? matches / reference.length : 0
  };
}

function needlemanSemiGlobal(reference, sequence, orientation) {
  const n = reference.length;
  const m = sequence.length;
  const width = m + 1;
  const total = (n + 1) * (m + 1);
  const score = new Int32Array(total);
  const trace = new Uint8Array(total);
  const gap = -3;
  const match = 2;
  const mismatch = -2;
  const idx = (i, j) => i * width + j;

  for (let i = 1; i <= n; i += 1) {
    score[idx(i, 0)] = 0;
    trace[idx(i, 0)] = 1;
  }
  for (let j = 1; j <= m; j += 1) {
    score[idx(0, j)] = 0;
    trace[idx(0, j)] = 2;
  }

  for (let i = 1; i <= n; i += 1) {
    for (let j = 1; j <= m; j += 1) {
      const diagonal = score[idx(i - 1, j - 1)] + (reference[i - 1] === sequence[j - 1] ? match : mismatch);
      const up = score[idx(i - 1, j)] + gap;
      const left = score[idx(i, j - 1)] + gap;
      let best = diagonal;
      let direction = 0;
      if (up > best) {
        best = up;
        direction = 1;
      }
      if (left > best) {
        best = left;
        direction = 2;
      }
      score[idx(i, j)] = best;
      trace[idx(i, j)] = direction;
    }
  }

  let bestI = n;
  let bestJ = m;
  let bestScore = -2147483648;
  for (let j = 0; j <= m; j += 1) {
    const value = score[idx(n, j)];
    if (value > bestScore) {
      bestScore = value;
      bestI = n;
      bestJ = j;
    }
  }
  for (let i = 0; i <= n; i += 1) {
    const value = score[idx(i, m)];
    if (value > bestScore) {
      bestScore = value;
      bestI = i;
      bestJ = m;
    }
  }

  const refOut = [];
  const readOut = [];
  let i = bestI;
  let j = bestJ;
  while (i > 0 && j > 0) {
    const direction = trace[idx(i, j)];
    if (direction === 0) {
      refOut.push(reference[i - 1]);
      readOut.push(sequence[j - 1]);
      i -= 1;
      j -= 1;
    } else if (direction === 1) {
      refOut.push(reference[i - 1]);
      readOut.push('-');
      i -= 1;
    } else {
      refOut.push('-');
      readOut.push(sequence[j - 1]);
      j -= 1;
    }
  }

  refOut.reverse();
  readOut.reverse();
  const alignedRef = refOut.join('');
  const alignedRead = readOut.join('');
  let matches = 0;
  let refBases = 0;
  for (let k = 0; k < alignedRef.length; k += 1) {
    if (alignedRef[k] !== '-') {
      refBases += 1;
      if (alignedRef[k] === alignedRead[k]) matches += 1;
    }
  }

  return {
    alignedRef,
    alignedRead,
    score: bestScore,
    orientation,
    refStartOffset: i,
    identity: refBases ? matches / refBases : 0
  };
}

function collectEdits(alignment) {
  const items = [];
  let refPos = alignment.refStartOffset || 0;
  let mismatchCount = 0;
  let insertionCount = 0;
  let deletionCount = 0;
  let coveredStart = 0;
  let coveredEnd = 0;
  let coveredStartColumn = -1;
  let coveredEndColumn = -1;

  for (let i = 0; i < alignment.alignedRef.length; i += 1) {
    const refBase = alignment.alignedRef[i];
    const readBase = alignment.alignedRead[i];
    if (refBase !== '-') {
      refPos += 1;
      if (readBase !== '-') {
        if (!coveredStart) coveredStart = refPos;
        coveredEnd = refPos;
        if (coveredStartColumn === -1) coveredStartColumn = i;
        coveredEndColumn = i;
      }
    }
  }

  refPos = alignment.refStartOffset || 0;

  for (let i = 0; i < alignment.alignedRef.length; i += 1) {
    const refBase = alignment.alignedRef[i];
    const readBase = alignment.alignedRead[i];
    if (refBase !== '-') refPos += 1;
    const insideCoveredReference = coveredStart && coveredEnd && refPos >= coveredStart && refPos <= coveredEnd;

    if (refBase !== '-' && readBase !== '-' && refBase !== readBase && insideCoveredReference) {
      mismatchCount += 1;
      items.push({ type: 'substitution', position: refPos, ref: refBase, alt: readBase });
    } else if (refBase !== '-' && readBase === '-' && insideCoveredReference) {
      deletionCount += 1;
      items.push({ type: 'deletion', position: refPos, ref: refBase, alt: '-' });
    } else if (refBase === '-' && readBase !== '-') {
      const position = Math.max(1, refPos);
      if (coveredStartColumn !== -1 && i > coveredStartColumn && i < coveredEndColumn) {
        insertionCount += 1;
        items.push({ type: 'insertion', position, ref: '+', alt: readBase });
      }
    }
  }

  const signature = items.slice(0, 12).map((edit) => {
    if (edit.type === 'substitution') return `${edit.ref}${edit.position}${edit.alt}`;
    if (edit.type === 'deletion') return `del${edit.position}`;
    return `ins${edit.position}${edit.alt}`;
  }).join('; ');
  const overflow = items.length > 12 ? `; +${items.length - 12} more` : '';

  return {
    items,
    mismatchCount,
    insertionCount,
    deletionCount,
    hasIndel: insertionCount + deletionCount > 0,
    hasAnyEdit: items.length > 0,
    signature: signature ? signature + overflow : '',
    coveredStart,
    coveredEnd
  };
}

function reverseComplement(sequence) {
  const map = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
  let out = '';
  for (let i = sequence.length - 1; i >= 0; i -= 1) out += map[sequence[i]] || 'N';
  return out;
}

function buildRunResult(samples, reference, settings, assay, config, inputFile) {
  const mutationRows = [];
  const alleleRows = [];
  const qcRows = [];
  for (const sample of samples) {
    for (const row of sample.positionRows) mutationRows.push({ ...row, label: sample.sample });
    for (const allele of sample.alleles) alleleRows.push({ ...allele, label: sample.sample });
    qcRows.push({ sample: sample.sample, label: sample.sample, ...sample.qc, files: sample.files.join('; '), warnings: sample.warnings.join('; ') });
  }
  const conversionRows = buildSubstitutionMatrix(samples);
  const targetRows = buildTargetRows(mutationRows, assay.targets, config.assay?.expectedEdit || null);
  return {
    createdAt: new Date().toISOString(),
    runName: config.runName || 'Amplicon analysis run',
    dataset: config.dataset || {},
    inputFile,
    reference,
    settings,
    assay: {
      targetSource: assay.source,
      targetWarning: assay.warning,
      targetPositions: Array.from(assay.targets).sort((a, b) => a - b),
      details: assay.details
    },
    samples,
    mutationRows,
    alleleRows,
    qcRows,
    conversionRows,
    targetRows,
    summary: summarize(samples, mutationRows, targetRows, conversionRows)
  };
}

function buildSubstitutionMatrix(samples) {
  const matrix = new Map();
  for (const sample of samples) {
    for (const row of sample.positionRows) {
      for (const alt of BASES) {
        const reads = row[alt] || 0;
        if (!reads || alt === row.refBase) continue;
        const key = `${sample.sample}|${row.refBase}|${alt}`;
        const current = matrix.get(key) || { sample: sample.sample, refBase: row.refBase, altBase: alt, reads: 0, denominator: 0 };
        current.reads += reads;
        current.denominator += row.denominator;
        matrix.set(key, current);
      }
    }
  }
  return Array.from(matrix.values()).map((row) => ({
    ...row,
    percent: row.denominator ? row.reads / row.denominator * 100 : 0
  })).sort((a, b) => b.reads - a.reads);
}

function buildTargetRows(mutationRows, targets, expectedEdit) {
  const positions = targets.size ? targets : new Set(mutationRows.map((row) => row.position));
  return mutationRows.filter((row) => positions.has(row.position)).map((row) => {
    const expectedReads = expectedEdit?.from === row.refBase && BASES.includes(expectedEdit.to)
      ? row[expectedEdit.to]
      : 0;
    return {
      ...row,
      expectedEditReads: expectedReads,
      expectedEditPercent: row.denominator ? expectedReads / row.denominator * 100 : 0
    };
  });
}

function summarize(samples, mutationRows, targetRows, conversionRows) {
  const totalPassed = samples.reduce((sum, sample) => sum + sample.qc.passedReads, 0);
  const maxEdit = Math.max(0, ...samples.map((sample) => sample.maxPercent));
  const topPosition = mutationRows.slice().sort((a, b) => b.percent - a.percent)[0] || null;
  const topTargetPosition = targetRows.slice().sort((a, b) => b.percent - a.percent)[0] || null;
  const topConversion = conversionRows[0] || null;
  return {
    sampleCount: samples.length,
    totalPassedReads: totalPassed,
    medianDepth: median(samples.map((sample) => sample.qc.passedReads)),
    maxEditPercent: maxEdit,
    topPosition,
    topTargetPosition,
    topConversion
  };
}

function writeOutputs(resultDir, run) {
  fs.writeFileSync(path.join(resultDir, 'summary.json'), JSON.stringify(toSerializableRun(run), null, 2) + '\n');
  fs.writeFileSync(path.join(resultDir, 'processed_mutation_table.csv'), buildMutationCsv(run));
  fs.writeFileSync(path.join(resultDir, 'allele_spectrum_table.csv'), buildAlleleCsv(run));
  fs.writeFileSync(path.join(resultDir, 'qc_read_counts.csv'), buildQcCsv(run));
  fs.writeFileSync(path.join(resultDir, 'substitution_matrix.csv'), buildConversionCsv(run));
  fs.writeFileSync(path.join(resultDir, 'target_window_table.csv'), buildTargetCsv(run));
  fs.writeFileSync(path.join(resultDir, 'amplicon_mutation_heatmap.svg'), buildSvgFigure(run));
  fs.writeFileSync(path.join(resultDir, 'run_report.md'), buildReport(run));
}

function toSerializableRun(run) {
  return {
    createdAt: run.createdAt,
    runName: run.runName,
    dataset: run.dataset,
    inputFile: run.inputFile,
    referenceLength: run.reference.length,
    settings: run.settings,
    assay: run.assay,
    summary: run.summary,
    qcRows: run.qcRows,
    topMutationRows: run.mutationRows.slice().sort((a, b) => b.percent - a.percent).slice(0, 20),
    targetRows: run.targetRows,
    substitutionMatrix: run.conversionRows,
    topAlleles: run.alleleRows.slice(0, 20)
  };
}

function buildMutationCsv(run) {
  const header = ['sample', 'label', 'position', 'ref_base', 'substitution_reads', 'deletion_reads', 'insertion_reads', 'signal_reads', 'denominator', 'percent', 'A', 'C', 'G', 'T', 'N'];
  const rows = run.mutationRows.map((row) => [
    row.sample,
    row.label,
    row.position,
    row.refBase,
    row.substitutionReads,
    row.deletionReads,
    row.insertionReads,
    row.displayReads,
    row.denominator,
    row.percent,
    row.A,
    row.C,
    row.G,
    row.T,
    row.N
  ]);
  return toCsv(header, rows);
}

function buildAlleleCsv(run) {
  const header = ['sample', 'label', 'sequence', 'n_reads', 'ratio', 'orientation', 'identity', 'edit_signature', 'mismatch_count', 'insertion_count', 'deletion_count', 'has_indel'];
  const rows = run.alleleRows.map((row) => [
    row.sample,
    row.label,
    row.sequence,
    row.nReads,
    row.ratio,
    row.orientation,
    row.identity,
    row.editSignature,
    row.mismatchCount,
    row.insertionCount,
    row.deletionCount,
    row.hasIndel
  ]);
  return toCsv(header, rows);
}

function buildQcCsv(run) {
  const header = ['sample', 'label', 'files', 'estimated_records', 'parsed_records', 'passed_reads', 'aligned_reads', 'low_identity_reads', 'no_indel_reads', 'edited_reads', 'unique_reads', 'mean_quality', 'mean_identity', 'dropped_short', 'dropped_quality', 'dropped_n', 'malformed_records', 'read_limit_reached', 'warnings'];
  const rows = run.qcRows.map((row) => [
    row.sample,
    row.label,
    row.files,
    row.estimatedRecords,
    row.parsedRecords,
    row.passedReads,
    row.alignedReads,
    row.lowIdentityReads,
    row.noIndelReads,
    row.editedReads,
    row.uniqueReads,
    row.meanQuality,
    row.meanIdentity,
    row.droppedShort,
    row.droppedQuality,
    row.droppedN,
    row.malformedRecords,
    row.readLimitReached,
    row.warnings
  ]);
  return toCsv(header, rows);
}

function buildConversionCsv(run) {
  const header = ['sample', 'ref_base', 'alt_base', 'reads', 'denominator', 'percent'];
  const rows = run.conversionRows.map((row) => [row.sample, row.refBase, row.altBase, row.reads, row.denominator, row.percent]);
  return toCsv(header, rows);
}

function buildTargetCsv(run) {
  const header = ['sample', 'position', 'ref_base', 'substitution_reads', 'deletion_reads', 'insertion_reads', 'signal_reads', 'denominator', 'percent', 'expected_edit_reads', 'expected_edit_percent', 'A', 'C', 'G', 'T', 'N'];
  const rows = run.targetRows.map((row) => [
    row.sample,
    row.position,
    row.refBase,
    row.substitutionReads,
    row.deletionReads,
    row.insertionReads,
    row.displayReads,
    row.denominator,
    row.percent,
    row.expectedEditReads,
    row.expectedEditPercent,
    row.A,
    row.C,
    row.G,
    row.T,
    row.N
  ]);
  return toCsv(header, rows);
}

function toCsv(header, rows) {
  const encode = (value) => `"${String(value == null ? '' : value).replace(/"/g, '""')}"`;
  return [header.map(encode).join(','), ...rows.map((row) => row.map(encode).join(','))].join('\n') + '\n';
}

function buildReport(run) {
  const topTargets = run.targetRows.slice().sort((a, b) => b.percent - a.percent).slice(0, 8);
  const expectedEdit = run.assay.details.expectedEdit;
  const totalAlignedReads = run.samples.reduce((sum, sample) => sum + sample.qc.alignedReads, 0);
  return [
    `# ${run.runName}`,
    '',
    `- Created: ${run.createdAt}`,
    `- Dataset: ${run.dataset.name || 'unspecified'}`,
    `- Source: ${run.dataset.sourceDocsUrl || 'unspecified'}`,
    `- FASTQ: ${run.dataset.fastqUrl || path.basename(run.inputFile)}`,
    `- Reference length: ${run.reference.length} bp`,
    `- Samples: ${run.summary.sampleCount}`,
    `- QC-passed reads: ${run.summary.totalPassedReads}`,
    `- Alignment-passed reads: ${totalAlignedReads}`,
    `- Edit signal shown in heatmap: ${getSignalLabel(run.settings.signalMode)}`,
    `- Target source: ${run.assay.targetSource}`,
    `- Target positions: ${run.assay.targetPositions.join(', ') || 'none'}`,
    `- Minimum alignment identity: ${formatNumber(run.settings.minIdentity * 100, 1)}%`,
    `- Max edit percentage: ${formatNumber(run.summary.maxEditPercent, 3)}%`,
    '',
    '## QC Summary',
    '',
    '| sample | QC-passed reads | alignment-passed reads | low-identity excluded | unique reads | mean Q | mean identity | max edit % | warnings |',
    '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |',
    ...run.samples.map((sample) => `| ${sample.sample} | ${sample.qc.passedReads} | ${sample.qc.alignedReads} | ${sample.qc.lowIdentityReads} | ${sample.qc.uniqueReads} | ${formatNumber(sample.qc.meanQuality, 2)} | ${formatNumber(sample.qc.meanIdentity * 100, 2)}% | ${formatNumber(sample.maxPercent, 3)}% | ${sample.warnings.join('; ') || 'OK'} |`),
    '',
    '## Target Window',
    '',
    '| position | ref | signal % | substitution reads | expected edit reads | expected edit % |',
    '| ---: | --- | ---: | ---: | ---: | ---: |',
    ...topTargets.map((row) => `| ${row.position} | ${row.refBase} | ${formatNumber(row.percent, 3)}% | ${row.substitutionReads} | ${row.expectedEditReads} | ${formatNumber(row.expectedEditPercent, 3)}% |`),
    '',
    '## Substitution Spectrum',
    '',
    '| conversion | reads | percent |',
    '| --- | ---: | ---: |',
    ...run.conversionRows.slice(0, 8).map((row) => `| ${row.refBase}>${row.altBase} | ${row.reads} | ${formatNumber(row.percent, 3)}% |`),
    '',
    '## Interpretation Notes',
    '',
    `- Expected edit model: ${expectedEdit ? `${expectedEdit.name || 'edit'} (${expectedEdit.from}>${expectedEdit.to})` : 'not specified'}.`,
    '- Reads were QC-filtered, unique-read aggregated, aligned to the amplicon in both orientations, and summarized at 1-based amplicon coordinates.',
    '- Edit percentages use position-level covered reads as the denominator; terminal no-coverage gaps are not counted as deletion edits.',
    '- This benchmark is intended as a reproducible validation fixture, not as a claim of equivalence to every CRISPResso2 output table.'
  ].join('\n') + '\n';
}

function buildSvgFigure(run) {
  const targets = new Set(run.assay.targetPositions);
  const cols = Array.from({ length: run.reference.length }, (_, idx) => ({ position: idx + 1, refBase: run.reference[idx] }));
  const samples = run.samples.slice();
  const compact = samples.length > 20;
  const cellW = compact ? (cols.length > 90 ? 7 : cols.length > 62 ? 9 : 12) : (cols.length > 90 ? 10 : 15);
  const rowH = compact ? 15 : 20;
  const left = compact ? 116 : 154;
  const right = 76;
  const top = 48;
  const bottom = 74;
  const width = Math.max(860, left + cols.length * cellW + right);
  const height = top + samples.length * rowH + bottom;
  const signalLabel = getSignalLabel(run.settings.signalMode);

  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-label="Amplicon mutation heatmap">`;
  svg += '<rect width="100%" height="100%" fill="#ffffff"/>';
  svg += '<style>.title{font:700 13px Arial,sans-serif;fill:#18212b}.axis{font:10px Arial,sans-serif;fill:#4c5b68}.row{font:11px Arial,sans-serif;fill:#18212b}.tiny{font:9px Arial,sans-serif;fill:#687786}.warn{fill:#b42318}</style>';
  svg += `<text x="14" y="20" class="title">${escapeHtml(run.runName)}</text>`;
  svg += `<text x="14" y="38" class="axis">${escapeHtml(signalLabel)} - generated ${escapeHtml(run.createdAt)}</text>`;

  cols.forEach((col, idx) => {
    const showTick = col.position === 1 || col.position === cols.length || col.position % 10 === 0 || targets.has(col.position);
    if (!showTick) return;
    const x = left + idx * cellW + cellW / 2;
    const y = targets.has(col.position) ? 38 : 34;
    svg += `<text x="${x}" y="${y}" text-anchor="middle" class="${targets.has(col.position) ? 'warn' : 'axis'}">${col.position}</text>`;
  });

  samples.forEach((sample, rowIndex) => {
    const y = top + rowIndex * rowH;
    const low = sample.qc.passedReads < run.settings.minDepth;
    svg += `<text x="${left - 8}" y="${y + rowH * 0.72}" text-anchor="end" class="row">${escapeHtml(clip(sample.sample, compact ? 18 : 27))}</text>`;
    sample.positionRows.forEach((row, colIndex) => {
      const x = left + colIndex * cellW;
      const target = targets.has(row.position);
      const fill = low ? '#e1e6ec' : heatColor(row.percent);
      svg += `<rect x="${x}" y="${y}" width="${Math.max(3, cellW - 1)}" height="${rowH - 1}" fill="${fill}" stroke="${target ? '#b42318' : '#ffffff'}" stroke-width="${target ? 1.5 : 0.6}"><title>${escapeHtml(sample.sample)} position ${row.position}: ${formatNumber(row.percent, 2)}%</title></rect>`;
    });
    svg += `<text x="${left + cols.length * cellW + 8}" y="${y + rowH * 0.72}" class="axis">${formatNumber(sample.maxPercent, 1)}%</text>`;
  });

  svg += buildLegend(left, height - 42);
  svg += '</svg>';
  return svg;
}

function heatColor(percent) {
  if (percent >= 50) return '#0f4c81';
  if (percent >= 25) return '#3a8ac6';
  if (percent >= 10) return '#8fc4ee';
  if (percent >= 2) return '#d7e9fb';
  return '#f8fbff';
}

function buildLegend(x, y) {
  const values = [0, 2, 10, 25, 50];
  let svg = `<g><text x="${x}" y="${y - 10}" class="axis">edit percentage</text>`;
  values.forEach((value, idx) => {
    svg += `<rect x="${x + idx * 50}" y="${y}" width="40" height="12" fill="${heatColor(value)}" stroke="#d8e1ea"/>`;
    svg += `<text x="${x + idx * 50}" y="${y + 29}" class="tiny">${value}%</text>`;
  });
  svg += `<text x="14" y="${y + 29}" class="tiny">Red outline: target position; gray row: below minimum reads</text>`;
  svg += '</g>';
  return svg;
}

function getSignalLabel(mode) {
  if (mode === 'substitution') return 'Substitution-only edit rate';
  if (mode === 'indel') return 'Indel-only edit rate';
  return 'Substitution + indel edit rate';
}

function escapeHtml(value) {
  return String(value == null ? '' : value).replace(/[&<>"']/g, (char) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;'
  }[char]));
}

function clip(value, length) {
  const text = String(value || '');
  return text.length > length ? `${text.slice(0, length - 1)}...` : text;
}

function formatNumber(value, digits = 1) {
  const n = Number(value);
  return Number.isFinite(n) ? n.toFixed(digits) : '0.0';
}

function printSummary(run, resultDir) {
  const topTarget = run.summary.topTargetPosition;
  const topConversion = run.summary.topConversion;
  console.log(`Run: ${run.runName}`);
  console.log(`Results: ${resultDir}`);
  console.log(`QC-passed reads: ${run.summary.totalPassedReads}`);
  console.log(`Alignment-passed reads: ${run.samples.reduce((sum, sample) => sum + sample.qc.alignedReads, 0)}`);
  console.log(`Low-identity reads excluded: ${run.samples.reduce((sum, sample) => sum + sample.qc.lowIdentityReads, 0)}`);
  console.log(`Alignment-passed unique reads: ${run.samples.reduce((sum, sample) => sum + sample.qc.uniqueReads, 0)}`);
  console.log(`Mean identity: ${formatNumber(run.samples[0].qc.meanIdentity * 100, 3)}%`);
  console.log(`Max edit: ${formatNumber(run.summary.maxEditPercent, 3)}%`);
  if (topTarget) console.log(`Top target position: ${topTarget.position} ${topTarget.refBase}, ${formatNumber(topTarget.percent, 3)}%`);
  if (topConversion) console.log(`Top conversion: ${topConversion.refBase}>${topConversion.altBase}, ${topConversion.reads} reads (${formatNumber(topConversion.percent, 3)}% of covered base opportunities)`);
}

try {
  main();
} catch (error) {
  console.error(error.stack || error.message || String(error));
  process.exitCode = 1;
}
