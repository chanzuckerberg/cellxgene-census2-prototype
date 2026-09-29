import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { copyFile, mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

async function buildFixture(t, change = () => {}) {
  const directory = await mkdtemp(join(tmpdir(), 'census-dashboard-links-'));
  t.after(() => rm(directory, { recursive: true, force: true }));
  for (const file of ['build.mjs', 'template.html', 'style.css', 'app.js', 'palettes.json']) {
    await copyFile(new URL(file, import.meta.url), join(directory, file));
  }
  const data = JSON.parse(await readFile(new URL('snapshot.json', import.meta.url), 'utf8'));
  change(data);
  await writeFile(join(directory, 'snapshot.json'), JSON.stringify(data));
  return spawnSync(process.execPath, [join(directory, 'build.mjs')], { encoding: 'utf8' });
}

test('the snapshot builds with canonical Notion page URLs', async (t) => {
  const result = await buildFixture(t);
  assert.equal(result.status, 0, result.stderr);
});

for (const kind of ['epic', 'task', 'source']) {
  test(`the build rejects ${kind} links without /p/`, async (t) => {
    const result = await buildFixture(t, (data) => {
      if (kind === 'source') data.sources.outline = data.sources.outline.replace('/p/', '/');
      else {
        const record = kind === 'epic' ? data.epics[0] : data.tasks[0];
        record.url = record.url.replace('/p/', '/');
      }
    });
    assert.notEqual(result.status, 0);
    assert.match(result.stderr, /Use the Notion page URL/);
  });
}
