(() => {
  'use strict';
  const data = JSON.parse(document.getElementById('project-data').textContent);
  const palettes = {
    current: {label:'Current', description:'The original epic colours use solid fills and white labels.', ink:'#ffffff', colors:Object.fromEntries(data.epics.map((epic) => [epic.id, epic.color]))},
    ...JSON.parse(document.getElementById('palette-data').textContent),
  };
  const $ = (id) => document.getElementById(id);
  const esc = (text) => String(text).replace(/[&<>"']/g, (c) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const link = (url, title, cls = '') => `<a class="${cls}" href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(title)} <span aria-hidden="true">↗</span></a>`;
  const epicLink = (epic, compact = false) => `<a class="epic-link${compact ? ' compact' : ''}" href="${esc(epic.url)}" target="_blank" rel="noopener noreferrer" style="--epic-color:var(--epic-${epic.id})"><i class="epic-swatch" aria-hidden="true"></i>${esc(compact ? `Phase ${epic.id}` : epic.title)} <span aria-hidden="true">↗</span></a>`;
  const epicBands = (epics) => epics.map((epic) => `<span class="epic-band" style="--epic-color:var(--epic-${epic.id})" title="${esc(epic.title)}"><span>Phase ${esc(epic.id)}</span></span>`).join('');
  const date = (value, options = {}) => new Intl.DateTimeFormat('en-GB', {day:'numeric', month:'short', timeZone:'UTC', ...options}).format(new Date(value));
  const statusClass = (value) => value.toLowerCase() === 'done' ? 'done' : value.toLowerCase() === 'in progress' ? 'active' : value.toLowerCase() === 'in review' ? 'review' : '';
  const badge = (status) => `<span class="badge ${statusClass(status)}">${esc(status)}</span>`;
  let selected = data.iterations.some((item) => `#iteration-${item.id}` === location.hash) ? location.hash.slice(11) : '1';
  const tasksFor = (id) => data.tasks.filter((task) => task.iteration === id);
  const hasPalette = (id) => Object.hasOwn(palettes,id);
  let savedPalette;
  try { savedPalette = localStorage.getItem('census2-palette'); } catch { /* Storage can be unavailable in embedded previews. */ }
  const requestedPalette = new URL(location.href).searchParams.get('palette');
  $('palette-options').innerHTML = Object.entries(palettes).map(([id,palette]) => `<button type="button" data-palette="${esc(id)}" aria-pressed="false"><span class="palette-preview" aria-hidden="true">${['4','5','6','9'].map((epic) => `<i style="background:${palette.colors[epic]}"></i>`).join('')}</span>${esc(palette.label)}</button>`).join('');

  function selectPalette(id, persist = true) {
    if (!hasPalette(id)) return;
    const palette = palettes[id];
    const style = document.documentElement.style;
    Object.entries(palette.colors).forEach(([epic,color]) => style.setProperty(`--epic-${epic}`,color));
    style.setProperty('--epic-ink',palette.ink);
    document.querySelectorAll('[data-palette]').forEach((button) => button.setAttribute('aria-pressed', String(button.dataset.palette === id)));
    $('palette-description').textContent = palette.description;
    if (persist) {
      try { localStorage.setItem('census2-palette',id); } catch { /* The selection still applies without storage. */ }
      try {
        const url = new URL(location.href);
        url.searchParams.set('palette',id);
        history.replaceState(null,'',url);
      } catch { /* The selection still applies without browser history. */ }
    }
  }
  selectPalette(hasPalette(requestedPalette) ? requestedPalette : hasPalette(savedPalette) ? savedPalette : 'current',false);

  const themes = {
    system: {label:'System', description:'The original green theme follows your device’s light or dark setting.'},
    slate: {label:'Slate', description:'Cool grey panels with blue accents.'},
    sand: {label:'Sand', description:'Warm ivory panels with terracotta accents.'},
    midnight: {label:'Midnight', description:'Dark navy panels with periwinkle accents.'},
  };
  const hasTheme = (id) => Object.hasOwn(themes,id);
  let savedTheme;
  try { savedTheme = localStorage.getItem('census2-theme'); } catch { /* Storage can be unavailable in embedded previews. */ }
  const requestedTheme = new URL(location.href).searchParams.get('theme');
  $('theme-options').innerHTML = Object.entries(themes).map(([id,theme]) => `<button type="button" data-theme-option="${id}" aria-pressed="false"><i class="theme-preview ${id}" aria-hidden="true"></i>${theme.label}</button>`).join('');

  function selectTheme(id, persist = true) {
    if (!hasTheme(id)) return;
    document.documentElement.dataset.theme = id;
    document.querySelectorAll('[data-theme-option]').forEach((button) => button.setAttribute('aria-pressed', String(button.dataset.themeOption === id)));
    $('theme-description').textContent = themes[id].description;
    if (persist) {
      try { localStorage.setItem('census2-theme',id); } catch { /* The selection still applies without storage. */ }
      try {
        const url = new URL(location.href);
        url.searchParams.set('theme',id);
        history.replaceState(null,'',url);
      } catch { /* The selection still applies without browser history. */ }
    }
  }
  selectTheme(hasTheme(requestedTheme) ? requestedTheme : hasTheme(savedTheme) ? savedTheme : 'system',false);

  $('source-links').innerHTML = link(data.sources.outline, 'Presentation outline') + link(data.sources.project, 'Project') + link(data.sources.tracker, 'Tracker');
  $('snapshot-date').textContent = `Snapshot · ${date(data.capturedAt, {year:'numeric'})}`;
  const age = Math.floor((Date.now() - Date.parse(data.capturedAt)) / 86400000);
  $('snapshot-age').textContent = age > 0 ? `Captured ${age} day${age === 1 ? '' : 's'} ago. Task status does not refresh automatically.` : 'Task status does not refresh automatically.';
  const done = data.tasks.filter((task) => task.status === 'Done').length;
  $('summary').innerHTML = [
    ['6','Iterations','plus work after the first release'],
    [data.epics.length,'Linked epics','from the project tracker'],
    [data.tasks.length,'Linked tasks',`${done} marked Done in the tracker`],
    ['30 Sep','First proposed gate','End-to-end skeleton'],
  ].map(([n,title,sub]) => `<div class="metric"><div><p class="metric-label">${esc(title)}</p><small>${esc(sub)}</small></div><strong>${esc(n)}</strong></div>`).join('');
  $('schedule-note').textContent = data.scheduleNote;
  $('mapping-note').textContent = data.mappingNote;
  $('epic-legend').innerHTML = data.epics.map((epic) => epicLink(epic)).join('');

  const rangeStart = Date.parse('2026-09-01T00:00:00Z');
  const rangeEnd = Date.parse('2027-02-01T00:00:00Z');
  const position = (value) => 100 * (Date.parse(`${value}T00:00:00Z`) - rangeStart) / (rangeEnd - rangeStart);
  const months = [['2026-09-01','Sep 2026'],['2026-10-01','Oct'],['2026-11-01','Nov'],['2026-12-01','Dec'],['2027-01-01','Jan 2027']];
  const grid = months.map(([value]) => `transparent ${position(value)}%, var(--line) ${position(value)}%, var(--line) calc(${position(value)}% + 1px), transparent calc(${position(value)}% + 1px)`).join(', ');
  let prior = 0;
  $('timeline').innerHTML = `<div class="axis"><span>Iteration / outline status</span><div class="months">${months.map(([value,label]) => `<span style="position:absolute;left:${position(value)}%">${label}</span>`).join('')}</div></div>` + data.iterations.map((item) => {
    const epics = data.epics.filter((epic) => epic.iterations.includes(item.id));
    let track;
    if (item.date) {
      const end = position(item.date);
      track = `<span class="gate-band" style="--start:${prior}%;--width:${end-prior}%"><span class="gate-label">${esc(date(item.date))}</span><span class="epic-bands">${epicBands(epics)}</span></span>`;
      prior = end;
    } else {
      track = `<span class="undated-group"><span class="undated">${item.id === '0' ? '<i class="dot done"></i> Design and measurement complete' : 'No proposed date'}</span><span class="undated-epics">${epicBands(epics)}</span></span>`;
    }
    return `<button type="button" class="iteration-row" data-iteration="${item.id}" style="min-height:${item.date ? 48 + epics.length * 19 : 90}px" aria-pressed="false" aria-controls="iteration-detail" aria-label="${esc(item.id === 'later' ? 'Later' : `Iteration ${item.id}`)}: ${esc(item.title)}${item.date ? `, proposed gate ${date(item.date,{year:'numeric'})}` : ''}; related phases ${epics.map((epic) => esc(epic.id)).join(', ')}"><span class="row-label"><span class="iteration-number">${item.id === 'later' ? '+' : `0${item.id}`}</span><span><span class="row-title">${esc(item.title)}</span><span class="row-meta"><i class="dot ${item.status === 'Done' ? 'done' : item.status === 'In progress' ? 'current' : ''}"></i>${esc(item.status)}</span></span></span><span class="track" style="background:linear-gradient(to right,${grid})">${track}</span></button>`;
  }).join('');

  function renderTasks() {
    const search = $('search').value.trim().toLowerCase();
    const status = $('status').value;
    const all = $('scope').value === 'all';
    const tasks = data.tasks.filter((task) => (all || task.iteration === selected) && (status === 'all' || task.status === status) && `${task.title} ${task.owner}`.toLowerCase().includes(search));
    $('tasks-heading').textContent = all ? 'Tasks across all iterations' : 'Tasks in this iteration';
    $('task-count').textContent = `${tasks.length} of ${all ? data.tasks.length : tasksFor(selected).length} tasks`;
    $('empty').hidden = tasks.length > 0;
    $('empty').textContent = tasksFor(selected).length === 0 && !all ? 'No child tasks are mapped to this iteration yet. Open the linked epics for their deliverables.' : 'No tasks match these filters.';
    $('tasks').innerHTML = tasks.map((task) => {
      const epic = data.epics.find((item) => item.id === task.epic);
      return `<tr><td>${link(task.url,task.title)}</td><td>${epicLink(epic,true)}</td><td><button class="table-iteration" data-iteration="${task.iteration}" type="button">${task.iteration === 'later' ? 'Later' : `Iteration ${task.iteration}`}</button></td><td>${esc(task.owner)}</td><td>${badge(task.status)}</td></tr>`;
    }).join('');
  }

  function selectIteration(id, setHash = true) {
    const item = data.iterations.find((iteration) => iteration.id === id);
    if (!item) return;
    selected = id;
    if (setHash) {
      try { history.replaceState(null, '', `#iteration-${id}`); } catch { /* The embedded preview can restrict browser history. */ }
    }
    document.querySelectorAll('.iteration-row').forEach((row) => row.setAttribute('aria-pressed', String(row.dataset.iteration === id)));
    const tasks = tasksFor(id);
    const completed = tasks.filter((task) => task.status === 'Done').length;
    const epics = data.epics.filter((epic) => epic.iterations.includes(id));
    $('iteration-detail').innerHTML = `<div class="detail-top"><p class="eyebrow">${id === 'later' ? 'DEFERRED WORK' : `ITERATION ${id}`} / PROPOSED MAPPING</p>${badge(item.status)}</div><h2 id="detail-title">${esc(item.title)}</h2><p class="delivery">${esc(item.delivery)}</p><div class="gate-box"><span class="gate-icon" aria-hidden="true">◇</span><div><p class="eyebrow">GATE SUMMARY${item.date ? ` / PROPOSED ${date(item.date,{year:'numeric'}).toUpperCase()}` : ''}</p><p>${esc(item.gate)}</p></div></div><div class="epic-heading"><h3>Related epics</h3><small>Each link opens Notion</small></div><div class="epic-list">${epics.map((epic) => epicLink(epic)).join('')}</div><div class="detail-notes"><h3>Source notes and open work</h3><ul>${item.notes.map((note) => `<li>${esc(note)}</li>`).join('')}</ul></div><div class="detail-bottom"><span>${tasks.length ? `${completed} of ${tasks.length} mapped tasks are marked Done.` : 'No child tasks are mapped yet.'}</span>${link(data.sources.outline,'Read the outline')}</div>`;
    renderTasks();
  }

  document.addEventListener('click', (event) => {
    const theme = event.target.closest('[data-theme-option]');
    if (theme) selectTheme(theme.dataset.themeOption);
    const palette = event.target.closest('[data-palette]');
    if (palette) selectPalette(palette.dataset.palette);
    const row = event.target.closest('[data-iteration]');
    if (row) {
      selectIteration(row.dataset.iteration);
      if (row.classList.contains('table-iteration')) $('iteration-detail').focus();
    }
  });
  $('search').addEventListener('input',renderTasks);
  $('status').addEventListener('change',renderTasks);
  $('scope').addEventListener('change',renderTasks);
  $('reset').addEventListener('click',() => { $('search').value=''; $('status').value='all'; $('scope').value='iteration'; renderTasks(); });
  window.addEventListener('hashchange',() => selectIteration(location.hash.slice(11),false));
  selectIteration(selected,false);
})();
