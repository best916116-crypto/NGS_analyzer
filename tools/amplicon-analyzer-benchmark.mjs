#!/usr/bin/env node

import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';

const BASES = ['A', 'C', 'G', 'T', 'N'];
const READ_LIMIT_ALL = Number.MAX_SAFE_INTEGER;

function main() {
  const configPath = path.resolve(process.argv[2] || 'benchmarks/public/crispresso2_base_editor/config.json');
  const configDir = path.dirname(configPath);
  const config = JSON.parse(fs.readFileSync(configPath, 'utf8').replace(/^\uFEFF/, ''));
  const resultDir = path.resolve(configDir, config.outputs?.resultDir || 'results');
  fs.mkdirSync(resultDir, { recursive: true });

  const setting = config.settingCsv ? parseSettingCsv(path.resolve(configDir, config.settingCsv)) : null;
  const reference = normalizeSeq(config.reference || setting?.reference);
  if (reference.length < 20) throw new Error('Reference amplicon sequence must be at least 20 bp.');

  const settings = normalizeSettings(config.settings);
  const targetPositions = config.targetPositions || setting?.targetPositions || '';
  const configAssay = config.assay && Object.keys(config.assay).length ? config.assay : {};
  const assayInput = Object.keys(configAssay).length ? configAssay : (setting?.assay || {});
  if (!assayInput.expectedEdit && setting?.expectedEdit) assayInput.expectedEdit = setting.expectedEdit;
  const assay = deriveAssayTargets(reference, assayInput, targetPositions);
  const groups = collectInputGroups(config, configDir, setting);
  if (!groups.length) throw new Error('No FASTQ input files were found.');
  const expectedCount = Number(config.dataset?.expectedSampleCount || setting?.expectedSampleCount || 0);
  if (expectedCount && groups.length !== expectedCount) {
    console.warn(`Warning: expected ${expectedCount} sample group(s), found ${groups.length}.`);
  }

  const samples = groups.map((group, index) => {
    const parsedGroup = parseInputGroup(group, settings);
    const sample = analyzeSample({
      sample: group.sample,
      files: group.files,
      roles: group.roles
    }, parsedGroup, reference, settings);
    console.log(`[${index + 1}/${groups.length}] ${group.sample}: ${sample.qc.passedReads} QC-passed, ${sample.qc.alignedReads} aligned`);
    return sample;
  });

  const inputFileLabel = config.dataset?.folder || config.dataset?.file || setting?.settingCsv || '';
  const run = buildRunResult(samples, reference, settings, assay, { ...config, assay: assayInput }, inputFileLabel);
  writeOutputs(resultDir, run);
  printSummary(run, resultDir);
}

function normalizeSettings(settings = {}) {
  return {
    minQ: Math.max(0, Number(settings.minQ ?? 20) || 0),
    minLen: Math.max(1, Number(settings.minLen ?? 30) || 1),
    readLimit: normalizeReadLimit(settings.readLimit),
    minIdentity: normalizeMinIdentity(settings.minIdentity ?? 0.8),
    dropN: settings.dropN !== false,
    minDepth: Math.max(0, Number(settings.minDepth ?? 0) || 0),
    signalMode: ['all', 'substitution', 'indel'].includes(settings.signalMode) ? settings.signalMode : 'all'
  };
}

function normalizeReadLimit(value) {
  const text = String(value ?? '').trim().toLowerCase().replace(/,/g, '');
  if (!text || text === '0' || text === 'all' || text === 'unlimited') return READ_LIMIT_ALL;
  const limit = Number(text);
  return Number.isFinite(limit) && limit > 0 ? Math.max(1, Math.floor(limit)) : READ_LIMIT_ALL;
}

function hasFiniteReadLimit(value) {
  return Number.isFinite(value) && value < READ_LIMIT_ALL / 2;
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
  return path.basename(filePath)
    .replace(/\.(fastq|fq|fastjoin|fastqjoin|fqjoin|join|txt)(\.gz)?$/i, '')
    .replace(/_S\d+_L\d{3}_R[12]_001$/i, '')
    .replace(/_S\d+_L\d{3}$/i, '')
    .replace(/_R[12]_001$/i, '')
    .replace(/[_\-.]R[12]$/i, '')
    .replace(/[_\-.][12]$/i, '')
    .replace(/\s+/g, '_') || 'sample';
}

function collectInputGroups(config, configDir, setting) {
  const dataset = config.dataset || {};
  let filePaths = [];
  if (dataset.folder) {
    const folder = path.resolve(configDir, dataset.folder);
    filePaths = fs.readdirSync(folder)
      .map((name) => path.join(folder, name))
      .filter((filePath) => fs.statSync(filePath).isFile() && isFastqLike(filePath));
  } else if (Array.isArray(dataset.files)) {
    filePaths = dataset.files.map((file) => path.resolve(configDir, file));
  } else if (dataset.file || config.fastq) {
    filePaths = [path.resolve(configDir, dataset.file || config.fastq)];
  } else if (setting?.files?.length && setting.baseDir) {
    filePaths = setting.files.map((file) => path.resolve(setting.baseDir, file));
  }

  const preferredKind = normalizePreferredInput(dataset.preferredInput || setting?.preferredInput || '');
  if (preferredKind === 'joined') filePaths = filePaths.filter(isJoinedFastq);
  if (preferredKind === 'raw') filePaths = filePaths.filter((filePath) => !isJoinedFastq(filePath));

  const start = dataset.sampleStart == null ? (setting?.sampleStart ?? null) : Number(dataset.sampleStart);
  const end = dataset.sampleEnd == null ? (setting?.sampleEnd ?? null) : Number(dataset.sampleEnd);
  const explicitSingleSampleName = dataset.sampleName && filePaths.length === 1 ? dataset.sampleName : '';
  const groups = new Map();
  for (const filePath of filePaths) {
    if (!fs.existsSync(filePath)) throw new Error(`FASTQ file not found: ${filePath}`);
    const sample = explicitSingleSampleName || sampleNameFromFile(filePath);
    const number = sampleNumberFromName(sample);
    if (Number.isFinite(start) && (number == null || number < start)) continue;
    if (Number.isFinite(end) && (number == null || number > end)) continue;
    if (!groups.has(sample)) groups.set(sample, { sample, files: [], paths: [], roles: new Set() });
    const group = groups.get(sample);
    group.files.push(path.basename(filePath));
    group.paths.push(filePath);
    group.roles.add(inferReadRole(filePath));
  }
  return Array.from(groups.values())
    .sort(compareGroups)
    .map((group) => ({ ...group, roles: Array.from(group.roles) }));
}

function isFastqLike(filePath) {
  return /\.(fastq|fq|fastjoin|fastqjoin|fqjoin|join|txt)(\.gz)?$/i.test(path.basename(filePath));
}

function isJoinedFastq(filePath) {
  return /\.(fastjoin|fastqjoin|fqjoin|join|txt)(\.gz)?$/i.test(path.basename(filePath));
}

function normalizePreferredInput(value) {
  const text = String(value || '').trim().toLowerCase();
  if (/join/.test(text)) return 'joined';
  if (/raw|fastq|fq|gz/.test(text)) return 'raw';
  return text;
}

function inferReadRole(filePath) {
  const name = path.basename(filePath);
  if (/(^|[_\-.])R?1([_\-.]|$)/i.test(name) || /_R1_001/i.test(name)) return 'R1';
  if (/(^|[_\-.])R?2([_\-.]|$)/i.test(name) || /_R2_001/i.test(name)) return 'R2';
  return 'single';
}

function sampleNumberFromName(value) {
  const match = String(value || '').match(/\d+/);
  return match ? Number(match[0]) : null;
}

function compareGroups(a, b) {
  const an = sampleNumberFromName(a.sample);
  const bn = sampleNumberFromName(b.sample);
  if (an != null && bn != null && an !== bn) return an - bn;
  if (an != null && bn == null) return -1;
  if (an == null && bn != null) return 1;
  return a.sample.localeCompare(b.sample, undefined, { numeric: true, sensitivity: 'base' });
}

function parseSettingCsv(settingCsv) {
  const rows = fs.readFileSync(settingCsv, 'utf8')
    .replace(/^\uFEFF/, '')
    .split(/\r?\n/)
    .map(parseCsvLine)
    .filter((row) => row.some(Boolean));
  if (rows.some((row) => row[0] === 'run_align_mutations')) return parseRowBasedSettingRows(settingCsv, rows);
  return parseModernSettingRows(settingCsv, rows);
}

function parseRowBasedSettingRows(settingCsv, rows) {
  const args = new Map();
  const opts = new Map();
  let runRow = null;
  for (const row of rows) {
    if (row[0] === 'arg') args.set(row[1], row[2] || '');
    if (row[0] === 'opts') opts.set(row[1], row[2] || '');
    if (row[0] === 'run_align_mutations') runRow = row;
  }
  if (!runRow) throw new Error(`Row-based setting CSV has no run_align_mutations row: ${settingCsv}`);
  const optsKey = runRow[1] || '';
  const referenceKey = runRow[2] || Array.from(args.keys()).find((key) => /^aseq/i.test(key));
  const siteKey = runRow[3] || Array.from(args.keys()).find((key) => /^site/i.test(key));
  const reference = normalizeSeq(args.get(referenceKey));
  const site = normalizeSeq(args.get(siteKey));
  const siteIndex = reference.indexOf(site);
  if (reference.length < 20 || !site || siteIndex < 0) throw new Error(`Could not map row-based setting CSV reference/site: ${settingCsv}`);
  const optionText = opts.get(optsKey) || '';
  const lengthMatch = optionText.match(/--user_region_length\s+(\d+)/i);
  const offsetMatch = optionText.match(/--user_region_beg_offset\s+(\d+)/i);
  const siteStart = siteIndex + 1;
  let targetStart = siteStart;
  let targetEnd = siteStart + site.length - 1;
  if (lengthMatch && offsetMatch) {
    const regionLength = Number(lengthMatch[1]);
    const beginOffset = Number(offsetMatch[1]);
    targetStart = Math.max(1, siteStart - beginOffset);
    targetEnd = Math.min(reference.length, targetStart + regionLength - 1);
  }
  return {
    format: 'row-based',
    settingCsv,
    baseDir: path.dirname(settingCsv),
    reference,
    site,
    targetPositions: `${targetStart}-${targetEnd}`,
    files: runRow.slice(4).map((value) => value.trim()).filter(isFastqLike)
  };
}

const SETTING_ALIASES = {
  siteName: ['site_name', 'site', 'site_id', 'name', 'amplicon_name'],
  reference: ['reference_amplicon_sequence', 'reference_amplicon', 'reference', 'ref_seq', 'refseq', 'amplicon_sequence', 'amplicon', 'aseq'],
  targetPositions: ['target_positions', 'targeting_window', 'target_window', 'targets', 'manual_targets', 'editing_window', 'user_region'],
  assayType: ['assay_type', 'assay', 'system', 'nuclease_type'],
  spacer: ['spacer_sequence', 'spacer', 'guide_sequence', 'guide', 'grna', 'sgrna', 'protospacer'],
  pam: ['pam_sequence', 'pam', 'pam_pattern'],
  spacerWindow: ['spacer_window', 'guide_window', 'cas_window', 'editing_window_in_spacer'],
  taleLeft: ['tale_left_sequence', 'left_tale_sequence', 'left_binding_site', 'tale_left', 'left_binding'],
  taleRight: ['tale_right_sequence', 'right_tale_sequence', 'right_binding_site', 'tale_right', 'right_binding'],
  taleSpacer: ['tale_spacer_sequence', 'tale_spacer', 'talen_spacer', 'known_spacer_sequence'],
  talePadding: ['tale_padding', 'window_padding', 'padding'],
  expectedEdit: ['expected_edit', 'base_edit', 'conversion'],
  sampleStart: ['sample_start', 'first_sample_number', 'first_sample', 'start_sample'],
  sampleEnd: ['sample_end', 'last_sample_number', 'last_sample', 'end_sample'],
  expectedSampleCount: ['expected_sample_count', 'sample_count', 'plate_size'],
  preferredInput: ['preferred_input', 'input_type']
};

function parseModernSettingRows(settingCsv, rows) {
  if (!rows.length) throw new Error(`Setting CSV is empty: ${settingCsv}`);
  const header = rows[0].map(normalizeSettingKey);
  if (!hasAnyHeader(new Set(header), SETTING_ALIASES.reference)) {
    throw new Error(`Setting CSV needs a reference_amplicon_sequence column: ${settingCsv}`);
  }
  const sites = rows.slice(1)
    .map((row, idx) => parseModernSettingRow(row, header, idx, settingCsv))
    .filter((site) => site && site.reference);
  if (!sites.length) throw new Error(`Setting CSV has no usable site rows: ${settingCsv}`);
  const first = sites[0];
  return {
    format: 'modern',
    settingCsv,
    baseDir: path.dirname(settingCsv),
    reference: first.reference,
    site: first.siteName,
    targetPositions: first.targetPositions,
    assay: first.assay,
    expectedEdit: first.expectedEdit,
    sampleStart: first.sampleStart,
    sampleEnd: first.sampleEnd,
    expectedSampleCount: first.expectedSampleCount,
    preferredInput: first.preferredInput,
    files: []
  };
}

function parseModernSettingRow(row, header, index, settingCsv) {
  const get = (key) => getSettingCell(row, header, SETTING_ALIASES[key]);
  const reference = normalizeSeq(get('reference'));
  if (!reference) return null;
  if (reference.length < 20) throw new Error(`Setting CSV row ${index + 2} reference is shorter than 20 bp: ${settingCsv}`);
  const rawSpacer = normalizeSeq(get('spacer'));
  const values = {
    type: normalizeAssayType(get('assayType')),
    spacer: rawSpacer,
    pam: normalizePattern(get('pam')),
    spacerWindow: get('spacerWindow').trim(),
    left: normalizeSeq(get('taleLeft')),
    right: normalizeSeq(get('taleRight')),
    taleSpacer: normalizeSeq(get('taleSpacer')),
    padding: get('talePadding').trim()
  };
  values.type = inferAssayType(values);
  if (values.type === 'cas') {
    values.spacer = rawSpacer || values.spacer;
    values.taleSpacer = '';
  } else if (values.type === 'tale') {
    values.taleSpacer = values.taleSpacer || rawSpacer;
    values.spacer = '';
  }
  const manualTargets = get('targetPositions').trim();
  const derived = deriveAssayTargets(reference, values, manualTargets);
  return {
    siteName: get('siteName').trim() || `Site ${index + 1}`,
    reference,
    targetPositions: formatTargetSet(derived.targets) || manualTargets,
    assay: values,
    expectedEdit: get('expectedEdit').trim().toUpperCase(),
    sampleStart: parseOptionalWholeNumber(get('sampleStart')),
    sampleEnd: parseOptionalWholeNumber(get('sampleEnd')),
    expectedSampleCount: parseOptionalWholeNumber(get('expectedSampleCount')),
    preferredInput: get('preferredInput').trim().toLowerCase()
  };
}

function normalizeSettingKey(value) {
  return String(value || '').trim().toLowerCase().replace(/[^a-z0-9]/g, '');
}

function hasAnyHeader(headerSet, aliases) {
  return aliases.some((alias) => headerSet.has(normalizeSettingKey(alias)));
}

function getSettingCell(row, header, aliases) {
  for (const alias of aliases) {
    const idx = header.indexOf(normalizeSettingKey(alias));
    if (idx !== -1) return unwrapSpreadsheetText(row[idx]);
  }
  return '';
}

function unwrapSpreadsheetText(value) {
  const text = String(value == null ? '' : value).trim();
  const match = text.match(/^=\s*"((?:[^"]|"")*)"\s*$/);
  return match ? match[1].replace(/""/g, '"').trim() : text;
}

function normalizeAssayType(value) {
  const key = String(value || '').trim().toLowerCase().replace(/[^a-z0-9]/g, '');
  if (!key) return '';
  if (['cas', 'crispr', 'crisprcas', 'cas9', 'cas12', 'baseeditor', 'primeeditor'].includes(key)) return 'cas';
  if (['tale', 'talen', 'taleen'].includes(key)) return 'tale';
  if (['custom', 'manual', 'window'].includes(key)) return 'custom';
  return '';
}

function inferAssayType(values) {
  if (values.type) return values.type;
  if (values.pam) return 'cas';
  if (values.left || values.right || values.taleSpacer || values.spacer) return 'tale';
  return 'custom';
}

function parseOptionalWholeNumber(value) {
  const text = String(value || '').trim();
  if (!text) return null;
  const number = Number(text);
  return Number.isFinite(number) && Number.isInteger(number) && number >= 0 ? number : null;
}

function formatTargetSet(targets) {
  const arr = Array.from(targets || []).sort((a, b) => a - b);
  if (!arr.length) return '';
  const ranges = [];
  let start = arr[0];
  let prev = arr[0];
  for (let i = 1; i <= arr.length; i += 1) {
    const value = arr[i];
    if (value === prev + 1) {
      prev = value;
      continue;
    }
    ranges.push(start === prev ? String(start) : `${start}-${prev}`);
    start = value;
    prev = value;
  }
  return ranges.join(',');
}

function parseCsvLine(line) {
  const cells = [];
  let cell = '';
  let quoted = false;
  for (let i = 0; i < line.length; i += 1) {
    const char = line[i];
    if (char === '"') {
      if (quoted && line[i + 1] === '"') {
        cell += '"';
        i += 1;
      } else {
        quoted = !quoted;
      }
    } else if (char === ',' && !quoted) {
      cells.push(cell.trim());
      cell = '';
    } else {
      cell += char;
    }
  }
  cells.push(cell.trim());
  return cells;
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
    details: {},
    spacerRegion: null
  });
  const assayType = normalizeAssayType(assay.type) || 'custom';

  if (assayType === 'cas') {
    const spacer = normalizeSeq(assay.spacer || '');
    const pamPattern = normalizePattern(assay.pam || '');
    const windowText = assay.spacerWindow || '4-8';
    const spacerWindow = parseTargets(windowText);
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
      source: `Cas ${match.orientation} spacer match - PAM ${pamPattern ? (match.pamMatched ? 'matched' : 'not confirmed') : 'not specified'} - spacer window ${windowText}`,
      warning: match.pamMatched ? '' : 'Cas mode: spacer was found, but the PAM pattern was not confirmed.',
      details: {
        spacer,
        pamPattern,
        spacerWindow: windowText,
        orientation: match.orientation,
        spacerStart: match.start,
        spacerEnd: match.end,
        expectedEdit: assay.expectedEdit || null
      },
      spacerRegion: null
    };
  }

  if (assayType === 'tale') {
    const left = normalizeSeq(assay.left || assay.taleLeft || '');
    const right = normalizeSeq(assay.right || assay.taleRight || '');
    const spacer = normalizeSeq(assay.taleSpacer || assay.spacer || '');
    const padding = Math.max(0, Number(assay.padding || assay.talePadding) || 0);
    if (left && right) {
      const match = findTalePair(reference, left, right);
      if (match) {
        return {
          targets: manual,
          source: `TALE/TALEN binding sites matched - spacer ${match.spacerLength} bp${padding ? ` - padding ${padding} bp` : ''}`,
          warning: '',
          details: { left, right, spacerLength: match.spacerLength, padding, spacerStart: match.start, spacerEnd: match.end, expectedEdit: assay.expectedEdit || null },
          spacerRegion: { start: match.start, end: match.end, label: 'Spacer' }
        };
      }
      return fallback('TALE/TALEN mode: left/right binding sites were not found in the reference.');
    }
    if (spacer) {
      const match = findSequenceEitherStrand(reference, spacer);
      if (match) {
        return {
          targets: manual,
          source: `TALE spacer-only ${match.orientation} match`,
          warning: 'TALE/TALEN mode: the spacer annotation was inferred from spacer-only matching. Left/right binding sites are recommended when available.',
          details: { spacer, orientation: match.orientation, spacerStart: match.start, spacerEnd: match.end, expectedEdit: assay.expectedEdit || null },
          spacerRegion: { start: match.start, end: match.end, label: 'Spacer' }
        };
      }
      return fallback('TALE/TALEN mode: spacer sequence was not found in the reference.');
    }
    return fallback('TALE/TALEN mode: no binding site or spacer sequence was provided.');
  }

  return fallback('');
}

function positionsFromRange(start, end, length) {
  const targets = new Set();
  const first = Math.max(1, Math.min(start, end));
  const last = Math.min(length, Math.max(start, end));
  for (let pos = first; pos <= last; pos += 1) targets.add(pos);
  return targets;
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

function findSequenceEitherStrand(reference, sequence) {
  const direct = reference.indexOf(sequence);
  if (direct !== -1) return { start: direct + 1, end: direct + sequence.length, orientation: 'forward' };
  const rc = reverseComplement(sequence);
  const reverse = reference.indexOf(rc);
  if (reverse !== -1) return { start: reverse + 1, end: reverse + rc.length, orientation: 'reverse' };
  return null;
}

function findTalePair(reference, left, right) {
  const leftVariants = uniqueSeqs([left, reverseComplement(left)]);
  const rightVariants = uniqueSeqs([right, reverseComplement(right)]);
  const candidates = [];
  leftVariants.forEach((leftSeq) => {
    rightVariants.forEach((rightSeq) => {
      const leftSites = allOccurrences(reference, leftSeq);
      const rightSites = allOccurrences(reference, rightSeq);
      leftSites.forEach((leftStart) => {
        rightSites.forEach((rightStart) => {
          const leftEnd = leftStart + leftSeq.length - 1;
          const rightEnd = rightStart + rightSeq.length - 1;
          if (leftEnd < rightStart) addTaleCandidate(candidates, leftEnd + 1, rightStart - 1, reference.length);
          else if (rightEnd < leftStart) addTaleCandidate(candidates, rightEnd + 1, leftStart - 1, reference.length);
        });
      });
    });
  });
  candidates.sort((a, b) => b.score - a.score || Math.abs(18 - a.spacerLength) - Math.abs(18 - b.spacerLength));
  return candidates[0] || null;
}

function uniqueSeqs(values) {
  return Array.from(new Set(values.filter(Boolean)));
}

function allOccurrences(reference, query) {
  const positions = [];
  let start = reference.indexOf(query);
  while (start !== -1) {
    positions.push(start + 1);
    start = reference.indexOf(query, start + 1);
  }
  return positions;
}

function addTaleCandidate(candidates, start, end, length) {
  if (end < start) return;
  const spacerLength = end - start + 1;
  const rangeScore = spacerLength >= 8 && spacerLength <= 30 ? 100 : 0;
  candidates.push({
    start: Math.max(1, start),
    end: Math.min(length, end),
    spacerLength,
    score: rangeScore + Math.max(0, 40 - Math.abs(18 - spacerLength))
  });
}

function emptyParsedGroup() {
  return {
    estimatedRecords: 0,
    totalRecords: 0,
    keptRecords: 0,
    malformedRecords: 0,
    droppedShort: 0,
    droppedQuality: 0,
    droppedN: 0,
    readLimitReached: false,
    unique: new Map(),
    meanQualityParts: [],
    totalQuality: 0,
    qualityReads: 0,
    basesKept: 0,
    readLengths: [],
    warnings: [],
    pairStats: {
      pairedReads: 0,
      joinedPairs: 0,
      unjoinedPairs: 0,
      unpairedReads: 0
    }
  };
}

function mergeParsedGroup(target, parsed) {
  target.estimatedRecords += parsed.estimatedRecords;
  target.totalRecords += parsed.totalRecords;
  target.keptRecords += parsed.keptRecords;
  target.malformedRecords += parsed.malformedRecords;
  target.droppedShort += parsed.droppedShort;
  target.droppedQuality += parsed.droppedQuality;
  target.droppedN += parsed.droppedN;
  target.readLimitReached = target.readLimitReached || parsed.readLimitReached;
  target.totalQuality += parsed.totalQuality || 0;
  target.qualityReads += parsed.qualityReads || 0;
  target.basesKept += parsed.basesKept || 0;
  if (parsed.keptRecords) target.meanQualityParts.push({ mean: parsed.meanQuality, n: parsed.keptRecords });
  for (const value of parsed.readLengths || []) {
    if (value) target.readLengths.push(value);
  }
  (parsed.warnings || []).forEach((warning) => target.warnings.push(warning));
  if (parsed.pairStats) {
    target.pairStats.pairedReads += parsed.pairStats.pairedReads || 0;
    target.pairStats.joinedPairs += parsed.pairStats.joinedPairs || 0;
    target.pairStats.unjoinedPairs += parsed.pairStats.unjoinedPairs || 0;
    target.pairStats.unpairedReads += parsed.pairStats.unpairedReads || 0;
  }
  parsed.unique.forEach((value, sequence) => {
    const current = target.unique.get(sequence) || { sequence, count: 0, qualitySum: 0 };
    current.count += value.count;
    current.qualitySum += value.qualitySum;
    target.unique.set(sequence, current);
  });
  target.uniqueReads = target.unique.size;
  target.meanQuality = groupMeanQuality(target);
  target.meanReadLength = target.keptRecords ? target.readLengths.reduce((sum, value) => sum + value, 0) / target.keptRecords : 0;
  target.medianReadLength = median(target.readLengths);
}

function groupMeanQuality(group) {
  if (group.qualityReads) return group.totalQuality / group.qualityReads;
  const total = group.meanQualityParts.reduce((sum, part) => sum + part.n, 0);
  if (!total) return 0;
  return group.meanQualityParts.reduce((sum, part) => sum + part.mean * part.n, 0) / total;
}

function parseFastqText(text, settings) {
  const parsed = emptyParsedGroup();
  const lines = String(text || '').replace(/\r\n/g, '\n').replace(/\r/g, '\n').split('\n');
  parsed.estimatedRecords = Math.floor(lines.length / 4);
  const limit = Math.min(settings.readLimit, parsed.estimatedRecords || settings.readLimit);

  for (let i = 0; i + 3 < lines.length && parsed.totalRecords < limit; i += 4) {
    const header = lines[i];
    const seqLine = lines[i + 1];
    const plus = lines[i + 2];
    const qualLine = lines[i + 3];
    if (!header || header[0] !== '@' || !plus || plus[0] !== '+') {
      parsed.malformedRecords += 1;
      continue;
    }
    addReadToParsed(parsed, seqLine, qualLine, settings);
  }
  parsed.readLimitReached = parsed.estimatedRecords > limit;
  return finalizeParsed(parsed);
}

function parseFastqRecords(text, limit) {
  const lines = String(text || '').replace(/\r\n/g, '\n').replace(/\r/g, '\n').split('\n');
  const estimatedRecords = Math.floor(lines.length / 4);
  const maxRecords = Math.min(Math.max(1, Number(limit) || estimatedRecords || 1), estimatedRecords || Number(limit) || 1);
  const records = [];
  let malformedRecords = 0;
  let totalRecords = 0;
  for (let i = 0; i + 3 < lines.length && totalRecords < maxRecords; i += 4) {
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
    records.push({ id: normalizeReadId(header), sequence: seq, quality: qual.slice(0, seq.length) });
  }
  return { estimatedRecords, totalRecords, malformedRecords, readLimitReached: estimatedRecords > maxRecords, records };
}

function addReadToParsed(parsed, seqLine, qualLine, settings) {
  parsed.totalRecords += 1;
  const seq = normalizeSeq(seqLine);
  const qual = String(qualLine || '').trim();
  if (!seq || qual.length < seq.length) {
    parsed.malformedRecords += 1;
    return false;
  }
  if (seq.length < settings.minLen) {
    parsed.droppedShort += 1;
    return false;
  }
  if (settings.dropN && seq.includes('N')) {
    parsed.droppedN += 1;
    return false;
  }
  const meanQ = meanQuality(qual.slice(0, seq.length));
  parsed.totalQuality += meanQ;
  parsed.qualityReads += 1;
  if (meanQ < settings.minQ) {
    parsed.droppedQuality += 1;
    return false;
  }
  parsed.keptRecords += 1;
  parsed.basesKept += seq.length;
  parsed.readLengths.push(seq.length);
  const current = parsed.unique.get(seq) || { sequence: seq, count: 0, qualitySum: 0 };
  current.count += 1;
  current.qualitySum += meanQ;
  parsed.unique.set(seq, current);
  return true;
}

function finalizeParsed(parsed) {
  parsed.uniqueReads = parsed.unique.size;
  parsed.meanQuality = parsed.qualityReads ? parsed.totalQuality / parsed.qualityReads : 0;
  parsed.meanReadLength = parsed.keptRecords ? parsed.basesKept / parsed.keptRecords : 0;
  parsed.medianReadLength = median(parsed.readLengths);
  return parsed;
}

function parseInputGroup(group, settings) {
  const parsedGroup = emptyParsedGroup();
  const joinedPaths = group.paths.filter(isJoinedFastq);
  const rawPaths = group.paths.filter((filePath) => !isJoinedFastq(filePath));
  if (joinedPaths.length) {
    for (const filePath of joinedPaths) mergeParsedGroup(parsedGroup, parseFastqText(readFastq(filePath), settings));
    if (rawPaths.length) parsedGroup.warnings.push('raw R1/R2 files skipped because joined FASTQ files were also provided');
    return finalizeParsed(parsedGroup);
  }

  const r1Paths = rawPaths.filter((filePath) => inferReadRole(filePath) === 'R1');
  const r2Paths = rawPaths.filter((filePath) => inferReadRole(filePath) === 'R2');
  const singlePaths = rawPaths.filter((filePath) => inferReadRole(filePath) === 'single');
  if (r1Paths.length && r2Paths.length) {
    const paired = parsePairedFastqFiles(r1Paths, r2Paths, settings);
    mergeParsedGroup(parsedGroup, paired);
    if (!paired.pairStats.joinedPairs && paired.pairStats.pairedReads) {
      parsedGroup.warnings.push('no paired-end overlaps found; analyzed R1/R2 independently');
      for (const filePath of rawPaths) mergeParsedGroup(parsedGroup, parseFastqText(readFastq(filePath), settings));
    }
    for (const filePath of singlePaths) mergeParsedGroup(parsedGroup, parseFastqText(readFastq(filePath), settings));
  } else {
    for (const filePath of rawPaths) mergeParsedGroup(parsedGroup, parseFastqText(readFastq(filePath), settings));
  }
  return finalizeParsed(parsedGroup);
}

function parsePairedFastqFiles(r1Paths, r2Paths, settings) {
  const parsed = emptyParsedGroup();
  const r1ById = new Map();
  const pairScanLimit = hasFiniteReadLimit(settings.readLimit)
    ? Math.max(settings.readLimit * 2, settings.readLimit + 1000)
    : READ_LIMIT_ALL;
  for (const filePath of r1Paths) {
    const records = parseFastqRecords(readFastq(filePath), pairScanLimit);
    parsed.estimatedRecords += records.estimatedRecords;
    parsed.malformedRecords += records.malformedRecords;
    parsed.readLimitReached = parsed.readLimitReached || records.readLimitReached;
    records.records.forEach((record) => {
      if (!r1ById.has(record.id)) r1ById.set(record.id, record);
    });
  }
  let pairAttempts = 0;
  for (const filePath of r2Paths) {
    const records = parseFastqRecords(readFastq(filePath), pairScanLimit);
    parsed.estimatedRecords += records.estimatedRecords;
    parsed.malformedRecords += records.malformedRecords;
    parsed.readLimitReached = parsed.readLimitReached || records.readLimitReached;
    for (const r2 of records.records) {
      if (hasFiniteReadLimit(settings.readLimit) && parsed.pairStats.joinedPairs >= settings.readLimit) {
        parsed.readLimitReached = true;
        break;
      }
      const r1 = r1ById.get(r2.id);
      if (!r1) {
        parsed.pairStats.unpairedReads += 1;
        continue;
      }
      pairAttempts += 1;
      const joined = joinPairedReads(r1, r2);
      if (!joined) {
        parsed.pairStats.unjoinedPairs += 1;
        continue;
      }
      parsed.pairStats.joinedPairs += 1;
      addReadToParsed(parsed, joined.sequence, joined.quality, settings);
    }
  }
  parsed.pairStats.pairedReads = pairAttempts;
  if (parsed.pairStats.unjoinedPairs) parsed.warnings.push(`${parsed.pairStats.unjoinedPairs} paired-end read(s) did not meet overlap criteria`);
  if (parsed.pairStats.unpairedReads) parsed.warnings.push(`${parsed.pairStats.unpairedReads} read(s) had no mate`);
  return finalizeParsed(parsed);
}

function normalizeReadId(header) {
  return String(header || '')
    .replace(/^@/, '')
    .split(/\s+/)[0]
    .replace(/\/[12]$/, '')
    .replace(/[_\-.]R[12]$/, '');
}

function joinPairedReads(r1, r2) {
  const r2Seq = reverseComplement(r2.sequence);
  const r2Qual = String(r2.quality || '').split('').reverse().join('');
  return mergeOverlappingReads(r1.sequence, r1.quality, r2Seq, r2Qual);
}

function mergeOverlappingReads(seq1, qual1, seq2, qual2) {
  const minOverlap = Math.min(12, seq1.length, seq2.length);
  let best = null;
  for (let offset = 0; offset <= seq1.length - minOverlap; offset += 1) {
    const start = Math.max(0, offset);
    const end = Math.min(seq1.length, offset + seq2.length);
    const overlap = end - start;
    if (overlap < minOverlap) continue;
    let matches = 0;
    let mismatches = 0;
    let compared = 0;
    for (let pos = start; pos < end; pos += 1) {
      const a = seq1[pos];
      const b = seq2[pos - offset];
      if (a === 'N' || b === 'N') continue;
      compared += 1;
      if (a === b) matches += 1;
      else mismatches += 1;
    }
    const identity = compared ? matches / compared : 0;
    const score = matches * 2 - mismatches * 5 + overlap * 0.03;
    if (identity >= 0.92 && (!best || score > best.score)) best = { offset, score, identity, overlap, mismatches };
  }
  if (!best) return null;
  return buildConsensus(seq1, qual1, seq2, qual2, best.offset);
}

function buildConsensus(seq1, qual1, seq2, qual2, offset) {
  const start = Math.min(0, offset);
  const end = Math.max(seq1.length, offset + seq2.length);
  let sequence = '';
  let quality = '';
  for (let pos = start; pos < end; pos += 1) {
    const i1 = pos;
    const i2 = pos - offset;
    const has1 = i1 >= 0 && i1 < seq1.length;
    const has2 = i2 >= 0 && i2 < seq2.length;
    if (has1 && has2) {
      const b1 = seq1[i1];
      const b2 = seq2[i2];
      const q1 = qScore(qual1[i1]);
      const q2 = qScore(qual2[i2]);
      if (b1 === b2) {
        sequence += b1;
        quality += qChar(Math.max(q1, q2));
      } else if (b1 === 'N') {
        sequence += b2;
        quality += qChar(q2);
      } else if (b2 === 'N') {
        sequence += b1;
        quality += qChar(q1);
      } else if (q2 >= q1) {
        sequence += b2;
        quality += qChar(q2);
      } else {
        sequence += b1;
        quality += qChar(q1);
      }
    } else if (has1) {
      sequence += seq1[i1];
      quality += qual1[i1] || '!';
    } else if (has2) {
      sequence += seq2[i2];
      quality += qual2[i2] || '!';
    }
  }
  return { sequence, quality };
}

function qScore(char) {
  return Math.max(0, String(char || '!').charCodeAt(0) - 33);
}

function qChar(score) {
  return String.fromCharCode(Math.max(0, Math.min(93, Math.round(score))) + 33);
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
    nonIndelSubstitutionReads: 0,
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
  let noIndelSubstitutionEditedReads = 0;
  let noIndelCMutationReads = 0;
  let noIndelGMutationReads = 0;
  let noIndelCgMutationReads = 0;
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
    let hasSubstitutionEdit = false;
    let hasCMutation = false;
    let hasGMutation = false;
    for (const edit of edits.items) {
      const stat = stats[edit.position - 1];
      if (!stat) continue;
      if (edit.type === 'substitution') {
        hasSubstitutionEdit = true;
        if (stat.refBase === 'C') hasCMutation = true;
        if (stat.refBase === 'G') hasGMutation = true;
        stat.substitutionReads += row.count;
        if (!edits.hasIndel) stat.nonIndelSubstitutionReads += row.count;
        const base = BASES.includes(edit.alt) ? edit.alt : 'N';
        stat[base] += row.count;
      } else if (edit.type === 'deletion') {
        stat.deletionReads += row.count;
      } else if (edit.type === 'insertion') {
        stat.insertionReads += row.count;
      }
      editedPositions.add(edit.position);
    }
    if (!edits.hasIndel) {
      if (hasSubstitutionEdit) noIndelSubstitutionEditedReads += row.count;
      if (hasCMutation) noIndelCMutationReads += row.count;
      if (hasGMutation) noIndelGMutationReads += row.count;
      if (hasCMutation || hasGMutation) noIndelCgMutationReads += row.count;
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
      nonIndelSubstitutionReads: stat.nonIndelSubstitutionReads,
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
  (parsedGroup.warnings || []).forEach((warning) => warnings.push(warning));
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
      noIndelSubstitutionEditedReads,
      noIndelCMutationReads,
      noIndelGMutationReads,
      noIndelCgMutationReads,
      uniqueReads: alleles.length,
      meanQuality: parsedGroup.meanQuality,
      medianReadLength: parsedGroup.medianReadLength,
      droppedShort: parsedGroup.droppedShort,
      droppedQuality: parsedGroup.droppedQuality,
      droppedN: parsedGroup.droppedN,
      malformedRecords: parsedGroup.malformedRecords,
      readLimitReached: parsedGroup.readLimitReached,
      pairedReads: parsedGroup.pairStats ? parsedGroup.pairStats.pairedReads : 0,
      joinedPairs: parsedGroup.pairStats ? parsedGroup.pairStats.joinedPairs : 0,
      unjoinedPairs: parsedGroup.pairStats ? parsedGroup.pairStats.unjoinedPairs : 0,
      unpairedReads: parsedGroup.pairStats ? parsedGroup.pairStats.unpairedReads : 0,
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
      details: assay.details,
      spacerRegion: assay.spacerRegion
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
  if (!targets.size) return [];
  return mutationRows.filter((row) => targets.has(row.position)).map((row) => {
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
  fs.writeFileSync(path.join(resultDir, 'summary.mutations.csv'), buildLegacyMutationSummaryCsv(run));
  fs.writeFileSync(path.join(resultDir, 'summary.nonX_per_pos.mutations.csv'), buildLegacyPerPositionMutationCsv(run));
  fs.writeFileSync(path.join(resultDir, 'amplicon_mutation_heatmap.svg'), buildSvgFigure(run));
  fs.writeFileSync(path.join(resultDir, 'run_report.md'), buildReport(run));
}

function toSerializableRun(run) {
  const summary = {
    ...run.summary,
    topPosition: run.summary.topPosition ? publicMutationRow(run.summary.topPosition) : null,
    topTargetPosition: run.summary.topTargetPosition ? publicMutationRow(run.summary.topTargetPosition) : null
  };
  return {
    createdAt: run.createdAt,
    runName: run.runName,
    dataset: run.dataset,
    inputFile: run.inputFile,
    referenceLength: run.reference.length,
    settings: run.settings,
    assay: run.assay,
    summary,
    qcRows: run.qcRows,
    topMutationRows: run.mutationRows.slice().sort((a, b) => b.percent - a.percent).slice(0, 20).map(publicMutationRow),
    targetRows: run.targetRows.map(publicMutationRow),
    substitutionMatrix: run.conversionRows,
    topAlleles: run.alleleRows.slice(0, 20)
  };
}

function publicMutationRow(row) {
  const { nonIndelSubstitutionReads, ...rest } = row;
  return rest;
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

function buildLegacyMutationSummaryCsv(run) {
  const header = ['filename', 'C', 'G', 'any', 'n_aligned_noindel', 'n_aligned_total'];
  const rows = run.samples.map((sample) => [
    legacyRowName(sample),
    legacyReadMetric(sample, 'noIndelCMutationReads', () => legacyReadMetric(sample, 'cMutationReads', () => fallbackReferenceMutationReads(sample, 'C'))),
    legacyReadMetric(sample, 'noIndelGMutationReads', () => legacyReadMetric(sample, 'gMutationReads', () => fallbackReferenceMutationReads(sample, 'G'))),
    legacyReadMetric(sample, 'noIndelCgMutationReads', () => legacyReadMetric(sample, 'cgMutationReads', () => legacyReadMetric(sample, 'noIndelSubstitutionEditedReads', () => sample.qc.editedReads || 0))),
    sample.qc.noIndelReads || 0,
    sample.qc.alignedReads || 0
  ]);
  return toCsv(header, rows);
}

function buildLegacyPerPositionMutationCsv(run) {
  const positionLabels = run.reference.split('').map((base, idx) => `${base || 'N'}${idx + 1}`);
  const header = ['nonX_per_pos', ...positionLabels, 'n_mut_any', 'n_aligned_noindel', 'n_aligned_total'];
  const rows = run.samples.map((sample) => {
    const countsByLabel = new Map();
    for (const row of sample.positionRows) {
      countsByLabel.set(`${row.refBase}${row.position}`, legacyPositionMutationReads(row));
    }
    return [
      legacyRowName(sample),
      ...positionLabels.map((label) => countsByLabel.get(label) || 0),
      legacyReadMetric(sample, 'noIndelSubstitutionEditedReads', () => legacyReadMetric(sample, 'substitutionEditedReads', () => sample.qc.editedReads || 0)),
      sample.qc.noIndelReads || 0,
      sample.qc.alignedReads || 0
    ];
  });
  return toCsv(header, rows);
}

function legacyRowName(sample) {
  const fileName = sample.files && sample.files.length ? sample.files[0] : sample.sample;
  return path.basename(String(fileName || sample.sample));
}

function legacyReadMetric(sample, key, fallback) {
  const value = Number(sample.qc && sample.qc[key]);
  if (Number.isFinite(value)) return value;
  return fallback ? fallback() : 0;
}

function fallbackReferenceMutationReads(sample, refBase) {
  return sample.positionRows
    .filter((row) => row.refBase === refBase)
    .reduce((sum, row) => sum + legacyPositionMutationReads(row), 0);
}

function legacyPositionMutationReads(row) {
  const value = Number(row.nonIndelSubstitutionReads);
  if (Number.isFinite(value)) return value;
  return nonReferenceSubstitutionReads(row);
}

function nonReferenceSubstitutionReads(row) {
  return BASES.reduce((sum, base) => {
    if (base === row.refBase) return sum;
    return sum + (Number(row[base]) || 0);
  }, 0);
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
  const header = ['sample', 'label', 'files', 'estimated_records', 'parsed_records', 'passed_reads', 'aligned_reads', 'low_identity_reads', 'no_indel_reads', 'edited_reads', 'unique_allele_sequences', 'mean_quality', 'mean_identity', 'dropped_short', 'dropped_quality', 'dropped_n', 'malformed_records', 'paired_reads', 'joined_pairs', 'unjoined_pairs', 'unpaired_reads', 'read_limit_reached', 'warnings'];
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
    row.pairedReads || 0,
    row.joinedPairs || 0,
    row.unjoinedPairs || 0,
    row.unpairedReads || 0,
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
    `- Read limit per sample: ${formatReadLimit(run.settings.readLimit)}`,
    `- Minimum alignment identity: ${formatNumber(run.settings.minIdentity * 100, 1)}%`,
    `- Max edit percentage: ${formatNumber(run.summary.maxEditPercent, 3)}%`,
    '',
    '## QC Summary',
    '',
    '| sample | QC-passed reads | alignment-passed reads | low-identity excluded | joined pairs | unjoined pairs | unique allele sequences | mean Q | mean identity | max per-base edit % | warnings |',
    '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |',
    ...run.samples.map((sample) => `| ${sample.sample} | ${sample.qc.passedReads} | ${sample.qc.alignedReads} | ${sample.qc.lowIdentityReads} | ${sample.qc.joinedPairs || 0} | ${sample.qc.unjoinedPairs || 0} | ${sample.qc.uniqueReads} | ${formatNumber(sample.qc.meanQuality, 2)} | ${formatNumber(sample.qc.meanIdentity * 100, 2)}% | ${formatNumber(sample.maxPercent, 3)}% | ${sample.warnings.join('; ') || 'OK'} |`),
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
    '- Raw R1/R2 FASTQ files were joined by high-confidence overlap consensus when possible, then reads were QC-filtered, unique-read aggregated, aligned to the amplicon in both orientations, and summarized at 1-based amplicon coordinates.',
    '- Edit percentages use position-level covered reads as the denominator; terminal no-coverage gaps are not counted as deletion edits.',
    '- This benchmark is intended as a reproducible validation fixture, not as a claim of equivalence to every CRISPResso2 output table.'
  ].join('\n') + '\n';
}

function buildSvgFigure(run) {
  const targets = new Set(run.assay.targetPositions);
  const targetRanges = rangesFromPositions(targets);
  const spacerRegion = normalizeRegion(run.assay.spacerRegion, run.reference.length);
  const cols = Array.from({ length: run.reference.length }, (_, idx) => ({ position: idx + 1, refBase: run.reference[idx] }));
  const samples = run.samples.slice();
  const compact = samples.length > 20;
  const cellW = compact ? (cols.length > 90 ? 7 : cols.length > 62 ? 9 : 12) : (cols.length > 90 ? 10 : 15);
  const rowH = compact ? 15 : 20;
  const left = compact ? 118 : 156;
  const right = 78;
  const top = 38;
  const bottom = spacerRegion ? 86 : 68;
  const sequenceH = compact ? 24 : 28;
  const width = Math.max(860, left + cols.length * cellW + right);
  const height = top + samples.length * rowH + sequenceH + bottom;

  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-label="Amplicon mutation heatmap">`;
  svg += '<rect width="100%" height="100%" fill="#ffffff"/>';
  svg += '<style>.title{font:700 italic 24px Arial,sans-serif;fill:#000}.axis{font:700 12px Arial,sans-serif;fill:#000}.row{font:700 13px Arial,sans-serif;fill:#000}.base{font:700 12px Arial,sans-serif;fill:#000}.tiny{font:700 10px Arial,sans-serif;fill:#000}.legend{font:700 14px Arial,sans-serif;fill:#000}.spacer{font:700 11px Arial,sans-serif;fill:#168b2b}</style>';
  svg += `<text x="${width / 2}" y="24" text-anchor="middle" class="title">${escapeHtml(run.runName)}</text>`;

  samples.forEach((sample, rowIndex) => {
    const y = top + rowIndex * rowH;
    const low = sample.qc.passedReads < run.settings.minDepth;
    svg += `<text x="${left - 8}" y="${y + rowH * 0.72}" text-anchor="end" class="row">${escapeHtml(clip(sample.sample, compact ? 18 : 27))}</text>`;
    sample.positionRows.forEach((row, colIndex) => {
      const x = left + colIndex * cellW;
      const fill = low ? '#e1e6ec' : heatColor(row.percent);
      svg += `<rect x="${x}" y="${y}" width="${cellW}" height="${rowH}" fill="${fill}" stroke="#000" stroke-width=".75" shape-rendering="crispEdges"><title>${escapeHtml(sample.sample)} position ${row.position}: ${formatNumber(row.percent, 2)}%</title></rect>`;
    });
    svg += `<text x="${left + cols.length * cellW + 8}" y="${y + rowH * 0.72}" class="axis">${formatNumber(sample.maxPercent, 1)}%</text>`;
  });

  targetRanges.forEach((range) => {
    const x = left + (range.start - 1) * cellW;
    const w = (range.end - range.start + 1) * cellW;
    svg += `<rect class="target-outline" x="${x}" y="${top}" width="${w}" height="${samples.length * rowH}" fill="none" stroke="#e00000" stroke-width="2.1" shape-rendering="crispEdges"><title>Target position ${range.start === range.end ? range.start : `${range.start}-${range.end}`}</title></rect>`;
  });

  const sequenceY = top + samples.length * rowH + (compact ? 17 : 20);
  const sequenceBoxY = sequenceY - (compact ? 13 : 15);
  svg += `<text x="${left - 8}" y="${sequenceY}" text-anchor="end" class="axis">Reference</text>`;
  cols.forEach((col, idx) => {
    const x = left + idx * cellW + cellW / 2;
    if (cellW >= 7) svg += `<text x="${x}" y="${sequenceY}" text-anchor="middle" class="base">${escapeHtml(col.refBase)}</text>`;
  });
  if (spacerRegion) {
    const x = left + (spacerRegion.start - 1) * cellW;
    const w = (spacerRegion.end - spacerRegion.start + 1) * cellW;
    svg += `<rect class="spacer-outline" x="${x}" y="${sequenceBoxY}" width="${w}" height="${compact ? 18 : 21}" fill="none" stroke="#18a538" stroke-width="2.2" shape-rendering="crispEdges"><title>${escapeHtml(spacerRegion.label || 'Spacer')} ${spacerRegion.start}-${spacerRegion.end}</title></rect>`;
  }

  svg += buildLegend(left, height - 64, run.settings.signalMode, !!spacerRegion);
  svg += '</svg>';
  return svg;
}

function heatColor(percent) {
  if (percent >= 50) return '#0046d8';
  if (percent >= 25) return '#0575e6';
  if (percent >= 10) return '#53aef5';
  if (percent >= 2) return '#d6ebff';
  return '#ffffff';
}

function buildLegend(x, y, signalMode, showSpacer) {
  const gradientId = `heatGradient${Math.floor(x)}${Math.floor(y)}`;
  const legendW = 260;
  let svg = `<defs><linearGradient id="${gradientId}" x1="0%" y1="0%" x2="100%" y2="0%"><stop offset="0%" stop-color="#ffffff"/><stop offset="42%" stop-color="#53aef5"/><stop offset="70%" stop-color="#0575e6"/><stop offset="100%" stop-color="#0046d8"/></linearGradient></defs>`;
  svg += `<g><rect x="${x}" y="${y}" width="${legendW}" height="12" fill="url(#${gradientId})" stroke="#000" stroke-width=".8"/>`;
  [0, 10, 25, 50].forEach((value) => {
    const tx = x + (value / 50) * legendW;
    svg += `<text x="${tx}" y="${y + 30}" text-anchor="middle" class="tiny">${value}</text>`;
  });
  svg += `<text x="${x + legendW / 2}" y="${y + 52}" text-anchor="middle" class="legend">${escapeHtml(getFigureLegendLabel(signalMode))}</text>`;
  if (showSpacer) {
    const lx = x + legendW + 42;
    svg += `<rect class="spacer-legend" x="${lx}" y="${y - 1}" width="34" height="14" fill="none" stroke="#18a538" stroke-width="2"/>`;
    svg += `<text x="${lx + 42}" y="${y + 11}" class="spacer">Spacer</text>`;
  }
  svg += '</g>';
  return svg;
}

function normalizeRegion(region, length) {
  if (!region) return null;
  const start = Math.max(1, Math.min(Number(region.start) || 0, Number(region.end) || 0));
  const end = Math.min(length, Math.max(Number(region.start) || 0, Number(region.end) || 0));
  if (!start || !end || end < start) return null;
  return { start, end, label: region.label || 'Spacer' };
}

function rangesFromPositions(positions) {
  const values = Array.from(positions || []).sort((a, b) => a - b);
  const ranges = [];
  if (!values.length) return ranges;
  let start = values[0];
  let end = values[0];
  for (let i = 1; i <= values.length; i += 1) {
    const value = values[i];
    if (value === end + 1) {
      end = value;
      continue;
    }
    ranges.push({ start, end });
    start = value;
    end = value;
  }
  return ranges;
}

function getFigureLegendLabel(mode) {
  if (mode === 'substitution') return 'Base substitution (%)';
  if (mode === 'indel') return 'Indel (%)';
  return 'Edit (%)';
}

function getSignalLabel(mode) {
  if (mode === 'substitution') return 'Base substitution (%)';
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

function formatReadLimit(value) {
  return hasFiniteReadLimit(value) ? `${Math.round(Number(value)).toLocaleString()} reads/sample` : 'all reads';
}

function printSummary(run, resultDir) {
  const topTarget = run.summary.topTargetPosition;
  const topConversion = run.summary.topConversion;
  console.log(`Run: ${run.runName}`);
  console.log(`Results: ${resultDir}`);
  console.log(`QC-passed reads: ${run.summary.totalPassedReads}`);
  console.log(`Alignment-passed reads: ${run.samples.reduce((sum, sample) => sum + sample.qc.alignedReads, 0)}`);
  console.log(`Low-identity reads excluded: ${run.samples.reduce((sum, sample) => sum + sample.qc.lowIdentityReads, 0)}`);
  console.log(`Unique allele sequences: ${run.samples.reduce((sum, sample) => sum + sample.qc.uniqueReads, 0)}`);
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
