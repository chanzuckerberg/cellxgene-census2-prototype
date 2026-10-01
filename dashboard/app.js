(() => {
  'use strict';
  const data = JSON.parse(document.getElementById('project-data').textContent);
  const palettes = {
    current: {label:'Current', description:'The original workstream colours use solid fills and white labels.', ink:'#ffffff', colors:Object.fromEntries(data.epics.map((epic) => [epic.id, epic.color]))},
    ...JSON.parse(document.getElementById('palette-data').textContent),
  };
  const $ = (id) => document.getElementById(id);
  const esc = (text) => String(text).replace(/[&<>"']/g, (c) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const link = (url, title, cls = '') => `<a class="${cls}" href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(title)} <span aria-hidden="true">↗</span></a>`;
  const epicLink = (epic, compact = false) => `<a class="epic-link${compact ? ' compact' : ''}" title="${esc(epic.title)}" href="${esc(epic.url)}" target="_blank" rel="noopener noreferrer" style="--epic-color:var(--epic-${epic.id})"><i class="epic-swatch" aria-hidden="true"></i>${esc(epic.workstream)} <span aria-hidden="true">↗</span></a>`;
  const statusIcons = {
    'Done':'<circle cx="10" cy="10" r="7"/><path d="m6.5 10 2.3 2.3 4.7-4.8"/>',
    'In Progress':'<circle cx="10" cy="10" r="7"/><path d="M10 3a7 7 0 0 1 0 14Z" fill="currentColor" stroke="none"/>',
    'Backlog':'<circle cx="10" cy="10" r="7"/>',
    'In Review':'<path d="M2 10s3-5 8-5 8 5 8 5-3 5-8 5-8-5-8-5Z"/><circle cx="10" cy="10" r="2"/>',
  };
  const statusIcon = (status) => `<svg class="status-icon" viewBox="0 0 20 20" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">${Object.hasOwn(statusIcons,status) ? statusIcons[status] : '<circle cx="10" cy="10" r="7"/><path d="M8 7a2 2 0 1 1 3 1.7c-1 .5-1 1.3-1 1.8M10 13h0"/>'}</svg>`;
  const epicBands = (epics) => epics.map((epic) => `<button type="button" class="epic-band" data-workstream="${epic.id}" aria-pressed="false" aria-controls="workstream-drawer" aria-haspopup="dialog" aria-label="${esc(epic.workstream)} workstream; tracker status: ${esc(epic.status)}. Show workstream details" style="--epic-color:var(--epic-${epic.id})" title="${esc(epic.workstream)} — Tracker status: ${esc(epic.status)}"><span class="workstream-name">${esc(epic.workstream)}</span>${statusIcon(epic.status)}</button>`).join('');
  const date = (value, options = {}) => new Intl.DateTimeFormat('en-GB', {day:'numeric', month:'short', timeZone:'UTC', ...options}).format(new Date(value));
  const statusClass = (value) => value.toLowerCase() === 'done' ? 'done' : value.toLowerCase() === 'in progress' ? 'active' : value.toLowerCase() === 'in review' ? 'review' : '';
  const badge = (status) => `<span class="badge ${statusClass(status)}">${esc(status)}</span>`;
  const defaultIteration = (data.iterations.find((item) => item.status === 'In progress') || data.iterations[0]).id;
  let selected = data.iterations.some((item) => `#iteration-${item.id}` === location.hash) ? location.hash.slice(11) : defaultIteration;
  let selectedEpic = null;
  let workstreamTrigger = null;
  const drawer = $('workstream-drawer');
  const narrowScreen = matchMedia('(max-width:1199px)');

  function openDrawer() {
    document.body.classList.add('drawer-open');
    if (!drawer.open) {
      if (narrowScreen.matches) drawer.showModal();
      else drawer.show();
    }
  }

  function closeDrawer() {
    if (!drawer.open) return;
    drawer.close();
    document.body.classList.remove('drawer-open');
    workstreamTrigger?.focus({preventScroll:true});
  }

  function dismissWorkstream() {
    selectIteration(selected);
  }

  $('close-workstream').addEventListener('click',dismissWorkstream);
  drawer.addEventListener('cancel',(event) => { event.preventDefault(); dismissWorkstream(); });
  document.addEventListener('keydown',(event) => {
    if (event.key === 'Escape' && drawer.open) { event.preventDefault(); dismissWorkstream(); }
  });
  drawer.addEventListener('click',(event) => {
    const bounds = drawer.getBoundingClientRect();
    if (narrowScreen.matches && event.target === drawer && (event.clientX < bounds.left || event.clientX > bounds.right || event.clientY < bounds.top || event.clientY > bounds.bottom)) dismissWorkstream();
  });
  narrowScreen.addEventListener('change',() => {
    if (drawer.open) { drawer.close(); openDrawer(); }
  });
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
  selectTheme(hasTheme(requestedTheme) ? requestedTheme : hasTheme(savedTheme) ? savedTheme : 'midnight',false);

  $('source-links').innerHTML = [
    ['design', 'Design doc'], ['twoPager', 'Meg’s 2-pager'], ['slides', 'Slides'],
    ['outline', 'Presentation outline'], ['project', 'Project'], ['tracker', 'Tracker'],
  ].map(([key, label]) => link(data.sources[key], label)).join('');
  $('snapshot-date').textContent = `Snapshot · ${date(data.capturedAt, {year:'numeric'})}`;
  const age = Math.floor((Date.now() - Date.parse(data.capturedAt)) / 86400000);
  $('snapshot-age').textContent = age > 0 ? `Captured ${age} day${age === 1 ? '' : 's'} ago. Tracker fields refresh every six hours.` : 'Tracker fields refresh every six hours.';
  const count = (items, status) => items.filter((item) => item.status.toLowerCase() === status).length;
  const dated = data.iterations.filter((item) => item.date);
  const nextGate = dated.find((item) => item.status !== 'Done') || dated[dated.length - 1];
  $('summary').innerHTML = [
    [data.iterations.filter((item) => item.id !== 'later').length,'Iterations',`${count(data.iterations,'done')} done, ${count(data.iterations,'in progress')} in progress, plus work after the first release`],
    [data.epics.length,'Linked epics',`${count(data.epics,'done')} done, ${count(data.epics,'in progress')} in progress in the tracker`],
    [data.tasks.length,'Linked tasks',`${count(data.tasks,'done')} done, ${count(data.tasks,'in review')} in review, ${count(data.tasks,'in progress')} in progress`],
    [date(nextGate.date),'Next gate',`${nextGate.title}, proposed`],
  ].map(([n,title,sub]) => `<div class="metric"><div><p class="metric-label">${esc(title)}</p><small>${esc(sub)}</small></div><strong>${esc(n)}</strong></div>`).join('');
  $('schedule-note').textContent = data.scheduleNote;
  $('mapping-note').textContent = data.mappingNote;
  $('epic-legend').innerHTML = data.epics.map((epic) => epicLink(epic)).join('');
  $('workstream-status-key').innerHTML = '<span class="workstream-status-label">Workstream status</span>' + [...new Set(data.epics.map((epic) => epic.status))].map((status) => `<span>${statusIcon(status)}${esc(status)}</span>`).join('');

  const rangeStart = Date.parse('2026-09-01T00:00:00Z');
  const rangeEnd = Date.parse('2027-02-01T00:00:00Z');
  const position = (value) => 100 * (Date.parse(`${value}T00:00:00Z`) - rangeStart) / (rangeEnd - rangeStart);
  const months = [['2026-09-01','Sep 2026'],['2026-10-01','Oct'],['2026-11-01','Nov'],['2026-12-01','Dec'],['2027-01-01','Jan 2027']];
  const grid = months.map(([value]) => `transparent ${position(value)}%, var(--line) ${position(value)}%, var(--line) calc(${position(value)}% + 1px), transparent calc(${position(value)}% + 1px)`).join(', ');
  let prior = 0;
  const demoLabel = (item) => item.id === '0' ? 'Review' : item.id === 'later' ? 'Future deliverables' : item.id === '5' ? 'Proposed release' : 'Proposed demo';
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
    return `<div class="iteration-row" data-iteration="${item.id}" style="min-height:${item.date ? 48 + epics.length * 35 : 90}px"><button type="button" class="row-label" data-iteration="${item.id}" aria-pressed="false" aria-controls="iteration-detail" aria-describedby="demo-${item.id}" aria-label="${esc(item.id === 'later' ? 'Later' : `Iteration ${item.id}`)}: ${esc(item.title)}${item.date ? `, proposed gate ${date(item.date,{year:'numeric'})}` : ''}"><span class="iteration-number">${item.id === 'later' ? '+' : `0${item.id}`}</span><span><span class="row-title">${esc(item.title)}</span><span class="row-meta ${item.status === 'Done' ? 'done' : item.status === 'In progress' ? 'current' : ''}"><i class="dot"></i>${esc(item.status)}</span><span class="row-deliverable" id="demo-${item.id}"><span class="deliverable-label">Delivers</span>${esc(item.delivery)}</span></span></button><div class="track" style="background:linear-gradient(to right,${grid})">${track}</div></div>`;
  }).join('');

  function renderTasks() {
    const search = $('search').value.trim().toLowerCase();
    const status = $('status').value;
    const all = $('scope').value === 'all';
    const scopedTasks = all ? data.tasks : selectedEpic !== null ? data.tasks.filter((task) => task.epic === selectedEpic) : tasksFor(selected);
    const tasks = scopedTasks.filter((task) => (status === 'all' || task.status === status) && `${task.title} ${task.owner}`.toLowerCase().includes(search));
    $('tasks-heading').textContent = all ? 'Tasks across all iterations' : selectedEpic !== null ? `Tasks in ${data.epics.find((epic) => epic.id === selectedEpic).workstream}` : 'Tasks in this iteration';
    $('task-count').textContent = `${tasks.length} of ${scopedTasks.length} tasks`;
    $('empty').hidden = tasks.length > 0;
    $('empty').textContent = scopedTasks.length === 0 && !all ? selectedEpic !== null ? 'No child tasks are linked to this workstream in the snapshot. Open the epic in Notion for its deliverables.' : 'No child tasks are mapped to this iteration yet. Open the linked epics for their deliverables.' : 'No tasks match these filters.';
    $('tasks').innerHTML = tasks.map((task) => {
      const epic = data.epics.find((item) => item.id === task.epic);
      return `<tr><td>${link(task.url,task.title)}</td><td>${epicLink(epic,true)}</td><td><button class="table-iteration" data-iteration="${task.iteration}" type="button">${task.iteration === 'later' ? 'Later' : `Iteration ${task.iteration}`}</button></td><td>${esc(task.owner)}</td><td>${badge(task.status)}</td></tr>`;
    }).join('');
  }

  function selectIteration(id, setHash = true) {
    const item = data.iterations.find((iteration) => iteration.id === id);
    if (!item) return;
    closeDrawer();
    selected = id;
    if (selectedEpic !== null) $('scope').value = 'selection';
    selectedEpic = null;
    $('scope').options[0].textContent = 'Selected iteration';
    if (setHash) {
      try { history.replaceState(null, '', `#iteration-${id}`); } catch { /* The embedded preview can restrict browser history. */ }
    }
    updateHighlights();
    const tasks = tasksFor(id);
    const completed = tasks.filter((task) => task.status === 'Done').length;
    const epics = data.epics.filter((epic) => epic.iterations.includes(id));
    $('iteration-detail').innerHTML = `<div class="detail-top"><p class="eyebrow">${id === 'later' ? 'DEFERRED WORK' : `ITERATION ${id}`} / PROPOSED MAPPING</p>${badge(item.status)}</div><h2 id="detail-title">${esc(item.title)}</h2><p class="delivery">${esc(item.delivery)}</p><div class="gate-box"><span class="gate-icon" aria-hidden="true">◇</span><div><p class="eyebrow">GATE SUMMARY${item.date ? ` / PROPOSED ${date(item.date,{year:'numeric'}).toUpperCase()}` : ''}</p><p>${esc(item.gate)}</p></div></div><div class="epic-heading"><h3>Related workstreams</h3><small>Each link opens Notion</small></div><div class="epic-list">${epics.map((epic) => epicLink(epic)).join('')}</div><div class="detail-notes"><h3>Source notes and open work</h3><ul>${item.notes.map((note) => `<li>${esc(note)}</li>`).join('')}</ul></div><div class="detail-bottom"><span>${tasks.length ? `${completed} of ${tasks.length} mapped tasks are marked Done.` : 'No child tasks are mapped yet.'}</span>${link(data.sources.outline,'Read the outline')}</div>`;
    $('iteration-detail').querySelector('.delivery').insertAdjacentHTML('afterend', `<section class="iteration-deliverables" aria-labelledby="deliverables-title"><h3 id="deliverables-title">${item.id === '0' ? 'Deliverables' : 'Planned deliverables'}</h3><ul>${item.deliverables.map((text) => `<li>${esc(text)}</li>`).join('')}</ul><p class="demo-summary"><span class="deliverable-label">${demoLabel(item)}</span>${esc(item.demo)}</p></section>`);
    renderTasks();
  }

  function updateHighlights() {
    document.querySelectorAll('.iteration-row').forEach((row) => {
      const active = selectedEpic === null && row.dataset.iteration === selected;
      row.classList.toggle('is-selected',active);
      row.querySelector('.row-label').setAttribute('aria-pressed',String(active));
    });
    document.querySelectorAll('[data-workstream]').forEach((button) => button.setAttribute('aria-pressed',String(button.dataset.workstream === selectedEpic)));
  }

  function selectEpic(id, setHash = true) {
    const epic = data.epics.find((item) => item.id === id);
    if (!epic) return;
    selectedEpic = id;
    if (!epic.iterations.includes(selected)) selected = epic.iterations[0];
    $('scope').value = 'selection';
    $('scope').options[0].textContent = 'Selected workstream';
    if (setHash) {
      try { history.replaceState(null,'',`#workstream-${id}`); } catch { /* The embedded preview can restrict browser history. */ }
    }
    updateHighlights();
    const tasks = data.tasks.filter((task) => task.epic === id);
    const completed = tasks.filter((task) => task.status === 'Done').length;
    const iterations = data.iterations.filter((item) => epic.iterations.includes(item.id));
    $('workstream-detail').innerHTML = `<div class="detail-top"><p class="eyebrow">WORKSTREAM / TRACKER SNAPSHOT</p>${badge(epic.status)}</div><h2 id="workstream-title">${esc(epic.workstream)}</h2><p class="tracker-reference">Tracker epic: ${esc(epic.title)}</p><p class="delivery">The status applies to the whole epic, across all mapped iterations.</p><div class="workstream-actions">${link(epic.url,'Open epic in Notion')}<button type="button" class="table-iteration" data-iteration="${selected}">Back to ${selected === 'later' ? 'Later' : `iteration ${selected}`}</button></div><div class="epic-heading"><h3>Related iterations</h3><small>Proposed mapping</small></div><div class="epic-list">${iterations.map((item) => `<button type="button" class="table-iteration" data-iteration="${item.id}">${item.id === 'later' ? 'Later' : `Iteration ${item.id}`}: ${esc(item.title)}</button>`).join('')}</div><div class="epic-heading"><h3>Linked tasks (${tasks.length})</h3><small>${completed} marked Done</small></div>${tasks.length ? `<ul class="drawer-tasks">${tasks.map((task) => `<li>${link(task.url,task.title)}<div class="drawer-task-meta">${badge(task.status)}<span>${esc(task.owner)}</span><span>${task.iteration === 'later' ? 'Later' : `Iteration ${esc(task.iteration)}`}</span></div></li>`).join('')}</ul>` : '<p class="workstream-note">No child tasks are linked in this snapshot. Open the epic in Notion for its deliverables.</p>'}<p class="workstream-note">The saved status and task counts do not establish completion of an iteration gate.</p>`;
    renderTasks();
    openDrawer();
    drawer.scrollTop = 0;
  }

  document.addEventListener('click', (event) => {
    const theme = event.target.closest('[data-theme-option]');
    if (theme) selectTheme(theme.dataset.themeOption);
    const palette = event.target.closest('[data-palette]');
    if (palette) selectPalette(palette.dataset.palette);
    const workstream = event.target.closest('[data-workstream]');
    if (workstream) {
      workstreamTrigger = workstream;
      selected = workstream.closest('.iteration-row').dataset.iteration;
      selectEpic(workstream.dataset.workstream);
      return;
    }
    const row = event.target.closest('[data-iteration]');
    if (row) {
      const fromDrawer = drawer.contains(row);
      selectIteration(row.dataset.iteration);
      if (fromDrawer) document.querySelector(`.row-label[data-iteration="${selected}"]`).focus({preventScroll:true});
      else if (row.classList.contains('table-iteration')) $('iteration-detail').focus();
    }
  });
  $('search').addEventListener('input',renderTasks);
  $('status').addEventListener('change',renderTasks);
  $('scope').addEventListener('change',renderTasks);
  function clearTaskFilters() {
    $('search').value = '';
    $('status').value = 'all';
    $('scope').value = 'selection';
  }
  $('reset').addEventListener('click',() => { clearTaskFilters(); renderTasks(); });
  $('dashboard-home').addEventListener('click',(event) => {
    if (event.button !== 0 || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    clearTaskFilters();
    selectIteration(defaultIteration,false);
    workstreamTrigger = null;
    document.querySelector('.appearance-settings').open = false;
    document.querySelector('.source-note').open = true;
    document.querySelectorAll('.timeline-scroll,.table-scroll').forEach((element) => element.scrollTo({left:0,top:0,behavior:'instant'}));
    try {
      const url = new URL(location.href);
      url.hash = '';
      history.replaceState(null,'',url);
    } catch { /* The view still resets without browser history. */ }
    $('dashboard-home').focus({preventScroll:true});
    window.scrollTo({left:0,top:0,behavior:'instant'});
  });
  function selectFromHash() {
    const workstream = location.hash.startsWith('#workstream-') && data.epics.find((epic) => epic.id === location.hash.slice(12));
    if (workstream) {
      selectIteration(workstream.iterations.includes(selected) ? selected : workstream.iterations[0],false);
      workstreamTrigger = document.querySelector(`.iteration-row[data-iteration="${selected}"] [data-workstream="${workstream.id}"]`);
      selectEpic(workstream.id,false);
    }
    else selectIteration(data.iterations.some((item) => `#iteration-${item.id}` === location.hash) ? location.hash.slice(11) : defaultIteration,false);
  }
  window.addEventListener('hashchange',selectFromHash);
  selectFromHash();
})();
