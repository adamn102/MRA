<!-- SPDX-License-Identifier: MIT -->
# Hosting the portal privately

The portal is a folder of static files. Hosting is "serve the folder over HTTPS", and
making it private is "put an authentication gate in front of the folder". Nothing about
the application changes.

## The principle that governs every choice

**The privacy guarantee is architectural, not a hosting feature.** Submitted sequences
never leave the browser because all computation is client-side. A login gate controls
*who can open the page*; it does not and cannot affect the sequences, which are never
transmitted. Two consequences:

1. Authentication protects the unpublished atlas and keeps the demo off the public web.
   It is access control, not data protection.
2. The host must not inject anything off-origin — analytics, CDN scripts, error
   reporters. `_headers` (shipped in this repo) enforces a Content-Security-Policy that
   makes any off-origin reference fail rather than silently work. **After deploying,
   verify:** open the browser network panel, submit a sequence, and confirm no request
   leaves the origin. This is the same check `scripts/serve_dev.py` is built around.

If collaborators will paste **real patient-derived** sequences, confirm the hosting
choice with your institution's IRB / IBC / data-governance office. Even with no data
transmitted, some institutions have a view on where the tool is hosted.

## Recommended: Cloudflare Pages + Cloudflare Access (free, Google or GitHub login)

Free for up to 50 users, your collaborators log in with their own Google/GitHub
accounts, and the gate runs at the edge before any file is served.

1. **Put the repo on GitHub** (private is fine; the repo is the source, not the gate).
2. **Cloudflare Pages → Create project → connect the repo.** Build command: none.
   Output directory: the folder containing `index.html` (the repo root, or `mra-portal/`
   if you keep the structure in this archive). Cloudflare serves `_headers`
   automatically.
3. **Add a custom domain** (Pages → Custom domains). HTTPS is automatic and is required:
   the service worker (offline) and the SHA-256 checksum verification both need a secure
   context.
4. **Turn on Cloudflare Access** (Zero Trust → Access → Applications → Add a
   self-hosted application):
   - Application domain: the portal's domain.
   - Identity provider: add Google (or GitHub) under Settings → Authentication.
   - Policy: Action **Allow**, Include → **Emails** and paste your collaborators'
     addresses (or Include → Emails ending in `@your-university.edu`).
   - Session duration: a day or a week, to taste.
5. **Test in an incognito window.** You should hit a Google login before the page loads,
   and be refused if you log in with an address not on the list.

To add or remove a collaborator later, edit the Access policy's email list. No redeploy.

Cloudflare logs authentication events (who logged in, when) with ~24 h retention on the
free tier. That is access metadata, never sequence data.

## Netlify (similar; password protection is a paid feature)

Netlify reads the same `_headers` file. Drag-and-drop the folder, or connect the repo
with no build command and the publish directory set to the folder with `index.html`.
Access control (password or Netlify Identity SSO) is on a paid tier; for free SSO
Cloudflare Access is the better fit, and you can even put Cloudflare Access in front of a
Netlify site.

## GitHub Pages — why it does NOT meet the "private" requirement

On Free, Pro and Team, a Pages site is **public even when its source repository is
private** — a private repo gives you a public site at an obscure URL, which is not access
control. Private Pages sites exist only on **GitHub Enterprise Cloud** (organization
accounts), and there every viewer must have a GitHub account with read access to the
repo — no external collaborators without a seat. Use GitHub for the source and the build
trigger, and a real auth gate (Cloudflare Access) for privacy.

GitHub Pages also does not serve custom response headers, so the `_headers` CSP cannot be
applied there; the portal ships a restrictive CSP `<meta>` fallback in each page, but
header-level CSP (as on Cloudflare/Netlify) is stronger.

## Institutional hosting behind SSO / VPN

The most defensible option for patient-derived data: ask your IT group for static hosting
behind your university's single sign-on (Shibboleth/SAML) or VPN. You copy the folder to
the served directory and translate the `_headers` rules into the server config:

```nginx
# nginx equivalent of _headers
location / {
    add_header Content-Security-Policy "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; worker-src 'self'; form-action 'none'; base-uri 'self'; frame-ancestors 'none'" always;
    add_header Referrer-Policy "no-referrer" always;
    add_header X-Content-Type-Options "nosniff" always;
}
location /data/ {
    # immutable release dirs
    location ~ /data/v[0-9.]+/ { add_header Cache-Control "public, max-age=31536000, immutable"; }
}
```

## Deploying updates

Whatever the host: rebuild assets if the model changed
(`scripts/build_web_assets.py ... --sdrd ...`), **run the full test suite including the
NL4-3 round-trip**, then publish. On Cloudflare/Netlify connected to the repo, a push
redeploys automatically. Cut a new atlas version as a new `data/vX.Y/` directory — never
overwrite a released one — so a citation keeps resolving (see `docs/VERSIONING.md`).

## A note on the version gate and auth

The `?v=` version machinery and the auth gate are independent. A collaborator following
`https://your-domain/?v=1.0` authenticates once, then the portal resolves the version
exactly as it does locally. Shareable links the UI generates still carry the explicit
version; they just now sit behind the login.
