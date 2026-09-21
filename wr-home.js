(() => {
  const tabs = [...document.querySelectorAll('[role="tab"]')];
  const description = document.getElementById('example-description');
  const descriptions = {
    compare: '<strong>The default is the last-event-assisted win ratio.</strong><p>This example compares exercise training with usual care in the included high-risk HF-ACTION subset, stratified by age. <a href="https://lmaowisc.github.io/WR/articles/WR_test_rec.html">See the full analysis ↗</a></p>',
    regression: '<strong>A treatment comparison adjusted for age.</strong><p>The proportional win-fractions model assumes a win ratio that is constant over follow-up. This example uses time-constant treatment and age covariates. <a href="https://lmaowisc.github.io/WR/articles/PW_reg.html">Read the model assumptions ↗</a></p>'
  };
  function selectTab(tab) {
    tabs.forEach(item => {
      const active = item === tab;
      item.setAttribute('aria-selected', String(active));
      item.tabIndex = active ? 0 : -1;
      document.getElementById(item.getAttribute('aria-controls')).hidden = !active;
    });
    description.innerHTML = descriptions[tab.dataset.example];
  }
  tabs.forEach((tab, index) => {
    tab.addEventListener('click', () => selectTab(tab));
    tab.addEventListener('keydown', event => {
      let next;
      if (event.key === 'ArrowRight') next = (index + 1) % tabs.length;
      if (event.key === 'ArrowLeft') next = (index + tabs.length - 1) % tabs.length;
      if (event.key === 'Home') next = 0;
      if (event.key === 'End') next = tabs.length - 1;
      if (next === undefined) return;
      event.preventDefault();
      selectTab(tabs[next]);
      tabs[next].focus();
    });
  });
  let toastTimer;
  function notify(message) {
    const toast = document.getElementById('copy-status');
    toast.textContent = message;
    toast.classList.add('visible');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => toast.classList.remove('visible'), 3200);
  }
  async function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      try { await navigator.clipboard.writeText(text); return; } catch (_) { /* Try the selection fallback below. */ }
    }
    const area = document.createElement('textarea');
    area.value = text;
    area.style.position = 'fixed';
    area.style.top = '-10000px';
    document.body.appendChild(area);
    area.select();
    const copied = document.execCommand('copy');
    area.remove();
    if (!copied) throw new Error('Clipboard unavailable');
  }
  document.querySelectorAll('[data-copy]').forEach(button => {
    button.addEventListener('click', async () => {
      const target = document.getElementById(button.dataset.copy);
      const text = target.innerText;
      try {
        await copyText(text);
        notify('Copied to clipboard');
      } catch (_) {
        notify('Select the code and copy it with your keyboard.');
      } finally {
        button.focus({ preventScroll: true });
      }
    });
  });
})();
