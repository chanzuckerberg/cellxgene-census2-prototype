import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { copyFile, mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

async function buildFixture(t, change = () => {}) {
  const directory = await mkdtemp(join(tmpdir(), 'census-dashboard-links-'));
  t.after(() => rm(directory, { recursive: true, force: true }));
  for (const file of ['build.mjs', 'template.html', 'style.css', 'app.js', 'palettes.json', 'favicon.svg']) {
    await copyFile(new URL(file, import.meta.url), join(directory, file));
  }
  const data = JSON.parse(await readFile(new URL('snapshot.json', import.meta.url), 'utf8'));
  change(data);
  await writeFile(join(directory, 'snapshot.json'), JSON.stringify(data));
  const result = spawnSync(process.execPath, [join(directory, 'build.mjs')], { encoding: 'utf8' });
  if (result.status === 0) result.html = await readFile(join(directory, 'index.html'), 'utf8');
  return result;
}

test('the snapshot builds with canonical Linear and Notion URLs', async (t) => {
  const result = await buildFixture(t);
  assert.equal(result.status, 0, result.stderr);
  const icon = result.html.match(/<link rel="icon"[^>]+href="data:image\/svg\+xml,([^"]+)"/);
  assert.ok(icon, 'The standalone page needs an embedded SVG favicon');
  assert.equal(decodeURIComponent(icon[1]), (await readFile(new URL('favicon.svg', import.meta.url), 'utf8')).trim());
  assert.doesNotMatch(result.html, /\/\*__[A-Z]+__\*\//);
});

for (const [key, invalid] of [
  ['design', 'https://app.notion.com/3d816a75378981f18f55c2bae3256a54'],
  ['twoPager', 'https://docs.google.com.evil.example/document/d/test/edit'],
  ['slides', 'javascript:alert(1)'],
  ['twoPager', 'https://docs.google.com/presentation/d/test/edit'],
  ['queryDiagram', 'https://scaling-dollop-pzyy6pk.pages.github.io.evil.example/census-query-sequence.html'],
]) {
  test(`the build rejects an invalid ${key} link: ${invalid}`, async (t) => {
    const result = await buildFixture(t, (data) => { data.sources[key] = invalid; });
    assert.notEqual(result.status, 0);
  });
}

for (const [field, invalid] of [['demo', ''], ['deliverables', []], ['deliverables', [' ']]]) {
  test(`the build rejects an iteration with invalid ${field}: ${JSON.stringify(invalid)}`, async (t) => {
    const result = await buildFixture(t, (data) => { data.iterations[0][field] = invalid; });
    assert.notEqual(result.status, 0);
  });
}

for (const [kind, message] of [['epic', /Use a Linear issue or project URL/], ['task', /Use a Linear issue or project URL/], ['tracker', /Use a Linear issue or project URL/], ['outline', /Use the Notion page URL/]]) {
  test(`the build rejects a malformed ${kind} link`, async (t) => {
    const result = await buildFixture(t, (data) => {
      if (kind === 'outline') data.sources.outline = data.sources.outline.replace('/p/', '/');
      else if (kind === 'tracker') data.sources.tracker = data.sources.tracker.replace('linear.app/', 'linear.app.evil.example/');
      else {
        const record = kind === 'epic' ? data.epics[0] : data.tasks[0];
        record.url = record.url.replace('/biohub/', '/other/');
      }
    });
    assert.notEqual(result.status, 0);
    assert.match(result.stderr, message);
  });
}
