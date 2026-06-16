#!/usr/bin/env node

import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const workDir = path.join(os.tmpdir(), `amplicon-analyzer-paired-validation-${Date.now()}`);
fs.mkdirSync(workDir, { recursive: true });

const reference = buildReference();
const targetPositions = '75-90';
const editPosition = 82;
const reads = buildReads(reference, editPosition, 240);

const rawDir = path.join(workDir, 'raw');
const joinedDir = path.join(workDir, 'joined');
fs.mkdirSync(rawDir, { recursive: true });
fs.mkdirSync(joinedDir, { recursive: true });

fs.writeFileSync(path.join(rawDir, 'synthetic_S1_L001_R1_001.fastq'), reads.r1);
fs.writeFileSync(path.join(rawDir, 'synthetic_S1_L001_R2_001.fastq'), reads.r2);
fs.writeFileSync(path.join(joinedDir, 'synthetic.fastqjoin'), reads.joined);

const rawSummary = runAnalyzer('raw-paired', rawDir, 'raw');
const joinedSummary = runAnalyzer('joined-reference', joinedDir, 'joined');

const rawSample = rawSummary.qcRows[0];
const joinedSample = joinedSummary.qcRows[0];
const rawTarget = rawSummary.summary.topTargetPosition || rawSummary.summary.topPosition;
const joinedTarget = joinedSummary.summary.topTargetPosition || joinedSummary.summary.topPosition;

const checks = [
  ['raw paired-end reads joined', rawSample.joinedPairs === reads.count],
  ['raw passed reads equal joined passed reads', rawSample.passedReads === joinedSample.passedReads],
  ['raw aligned reads equal joined aligned reads', rawSample.alignedReads === joinedSample.alignedReads],
  ['expected edit position retained', rawTarget.position === editPosition && joinedTarget.position === editPosition],
  ['edit percentage matches joined FASTQ', nearlyEqual(rawTarget.percent, joinedTarget.percent, 1e-9)]
];

const failed = checks.filter(([, ok]) => !ok);
console.log(`Paired-end validation workspace: ${workDir}`);
checks.forEach(([name, ok]) => console.log(`${ok ? 'PASS' : 'FAIL'} ${name}`));
console.log(`Raw joined pairs: ${rawSample.joinedPairs}`);
console.log(`Raw aligned reads: ${rawSample.alignedReads}`);
console.log(`Joined aligned reads: ${joinedSample.alignedReads}`);
console.log(`Edit position ${editPosition}: raw ${rawTarget.percent.toFixed(6)}%, joined ${joinedTarget.percent.toFixed(6)}%`);

if (failed.length) process.exit(1);

function runAnalyzer(name, folder, preferredInput) {
  const resultDir = path.join(workDir, `${name}-results`);
  const configPath = path.join(workDir, `${name}.json`);
  fs.writeFileSync(configPath, JSON.stringify({
    runName: `paired-end validation ${name}`,
    reference,
    targetPositions,
    dataset: {
      folder,
      preferredInput,
      expectedSampleCount: 1
    },
    settings: {
      readLimit: reads.count,
      minQ: 20,
      minLen: 30,
      minIdentity: 0.8,
      dropN: true,
      signalMode: 'substitution'
    },
    outputs: { resultDir }
  }, null, 2));
  const result = spawnSync(process.execPath, [path.join(repoRoot, 'tools', 'amplicon-analyzer-benchmark.mjs'), configPath], {
    cwd: repoRoot,
    encoding: 'utf8'
  });
  if (result.status !== 0) {
    process.stdout.write(result.stdout || '');
    process.stderr.write(result.stderr || '');
    throw new Error(`Analyzer failed for ${name}`);
  }
  return JSON.parse(fs.readFileSync(path.join(resultDir, 'summary.json'), 'utf8'));
}

function buildReference() {
  const seed = 'ACGTTGCAAGTCGATCGTACCGATGCTAGTCAGTACGATCGATGCTACGATCGTAGCTAGCATCGATCGTACGATCGATGCTAGCTACGATCGTAGCTAGCATGCTAGTCGATCGTAGCTAGCATCGATGCTAGTCGATCGTAGCTAGCATCGATCGTACGATCGTAGCTAGCATGCTAGTCGATCGTAGCATCGATCGTACGATCGTAGCTA';
  return seed.slice(0, 220);
}

function buildReads(ref, editPos, count) {
  const r1 = [];
  const r2 = [];
  const joined = [];
  for (let i = 0; i < count; i += 1) {
    let seq = ref;
    if (i < count / 2) seq = replaceAt(seq, editPos - 1, ref[editPos - 1] === 'C' ? 'T' : 'C');
    const id = `@SYNTH:${i}:1:1`;
    const read1 = seq.slice(0, 150);
    const read2 = reverseComplement(seq.slice(seq.length - 150));
    r1.push(`${id} 1:N:0:1`, read1, '+', 'I'.repeat(read1.length));
    r2.push(`${id} 2:N:0:1`, read2, '+', 'I'.repeat(read2.length));
    joined.push(`${id} 1:N:0:1`, seq, '+', 'I'.repeat(seq.length));
  }
  return {
    count,
    r1: `${r1.join('\n')}\n`,
    r2: `${r2.join('\n')}\n`,
    joined: `${joined.join('\n')}\n`
  };
}

function reverseComplement(sequence) {
  const map = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
  return String(sequence).split('').reverse().map((base) => map[base] || 'N').join('');
}

function replaceAt(value, index, char) {
  return value.slice(0, index) + char + value.slice(index + 1);
}

function nearlyEqual(a, b, tolerance) {
  return Math.abs(Number(a) - Number(b)) <= tolerance;
}
