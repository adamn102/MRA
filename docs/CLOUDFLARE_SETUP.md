<!-- SPDX-License-Identifier: MIT -->
# Cloudflare Pages + Access (Google login) — setup runbook

End state: the portal served over HTTPS at a domain you choose, with a Google login in
front, open only to the email addresses you list. Free for up to 50 users.

Some of this I cannot do for you: it happens in *your* Cloudflare and Google accounts and
needs your credentials. Those steps are marked **[you, in the dashboard]**. The repo-side
pieces are already done and committed.

There are two independent halves: **(A) host the files** and **(B) put Google login in
front**. Do A first, confirm the site loads, then do B.

--------------------------------------------------------------------------------
## Part A — host the static files on Cloudflare Pages

### A1. Put the folder in a Git repo  [you]
Cloudflare Pages deploys from GitHub or GitLab. The repo may be **private** — Pages reads
it through an authorized connection, and the published site's privacy is handled by Access
in Part B, not by repo visibility.

Make the portal folder the thing Pages serves. Either:
- make the repo root the portal folder (so `index.html` is at the top), or
- keep this archive's layout and set the output directory to `mra-portal` in A3.

`_headers`, `_redirects` (none needed), `wrangler.toml` and the `data/` tree all travel
with the repo. Do **not** add a `.nojekyll`-style build; there is no build.

### A2. Create the Pages project  [you, in the dashboard]
Cloudflare dashboard → **Workers & Pages** → **Create** → **Pages** → **Connect to Git**
→ authorize and pick the repo.

### A3. Build settings  [you, in the dashboard]
- **Framework preset:** None.
- **Build command:** leave empty.
- **Build output directory:** the folder that contains `index.html` — `/` if the repo
  root is the portal, or `mra-portal` if you kept this archive's layout.
- Save and deploy. The first build just copies the folder; it finishes in seconds.

You now have a URL like `mra-portal.pages.dev`. **It is public at this point** — anyone
with the link can open it. That is expected; Part B closes it. Until Part B is done, do
not share the `.pages.dev` link or paste real sequences.

### A4. Custom domain (recommended, and needed for a clean Access setup)  [you]
Pages project → **Custom domains** → add e.g. `mra.yourlab.org`. HTTPS is automatic.
A custom domain is worth it here: the service worker (offline use) and the SHA-256
checksum verification both require a secure context, and Access policies are cleaner on a
domain you control. If your lab's DNS is not on Cloudflare, Cloudflare will show the CNAME
to add; or use the `.pages.dev` domain for Access directly, which also works.

### A5. Confirm the headers are live  [you]
Load the site, open the browser dev tools → Network → click the document request →
Response Headers. You should see `Content-Security-Policy` and `Referrer-Policy` exactly
as in `_headers`. Then the privacy check that matters: with the Network panel open, go to
the sequence tool, paste a sequence, run it, and confirm **no request leaves the origin**.

--------------------------------------------------------------------------------
## Part B — Google login in front, with Cloudflare Access

Access is part of Cloudflare **Zero Trust**. The free plan covers up to 50 users.

### B1. Open Zero Trust and pick a team name  [you, in the dashboard]
Dashboard → **Zero Trust**. First time only, it asks for a team name (becomes
`yourteam.cloudflareaccess.com`, the login domain your collaborators briefly see) and a
plan — choose **Free**. A card may be requested at onboarding; the free plan does not
charge for Access up to 50 users.

### B2. Add Google as an identity provider  [you, in the dashboard]
Zero Trust → **Settings** → **Authentication** → **Login methods** → **Add new**.

Two choices:
- **Google** (personal Gmail / any Google account) or **Google Workspace** (a whole
  `@yourlab.org` domain). Both require a Google Cloud OAuth client: in Google Cloud
  Console create an **OAuth 2.0 Client ID** (type: Web application), set the authorized
  redirect URI to the value Cloudflare shows on the Google login-method form
  (`https://<yourteam>.cloudflareaccess.com/cdn-cgi/access/callback`), then paste the
  client ID and secret back into Cloudflare. Cloudflare's form has the exact redirect URI
  to copy.
- **Zero-setup alternative: One-time PIN.** If creating a Google OAuth client is more
  than you want right now, enable **One-time PIN** instead (no IdP, no Google Cloud step):
  Access emails each collaborator a six-digit code. You still control exactly who gets in
  via the email list in B4. You can switch to Google later without disrupting anything.

Use **Test** on the login method to confirm it works before continuing.

### B3. Create the Access application  [you, in the dashboard]
Zero Trust → **Access** → **Applications** → **Add an application** → **Self-hosted**.
- **Application name:** anything neutral (e.g. "MRA portal").
- **Session duration:** a day or a week, to taste.
- **Public hostname:** the domain from A4 (`mra.yourlab.org`, or the `.pages.dev` host).
- **Identity providers:** tick Google (and/or One-time PIN).

### B4. Write the allow policy — email whitelist  [you, in the dashboard]
This is the whitelist. Add one policy on the application:
- **Action:** Allow.
- **Include** → **Emails** → paste your collaborators' addresses, one per line.
- **Require** → **Login methods** → **Google** (so the list is satisfied only by a
  Google login, not by one-time PIN or any other method you might add later).
- Leave everything else default. Save.

**How the two conditions combine:** a visitor gets in only if they (1) authenticate with
Google **and** (2) sign in with an address on the list. A Google account whose address is
not listed is refused; a listed address that cannot pass Google login is refused. That is
exactly "only whitelisted Google addresses."

Add a second policy below it — **Action: Block, Include: Everyone** — so the default is
deny. (Access denies anything not explicitly allowed anyway, but an explicit block policy
makes the intent unmistakable to the next person who reads it.)

### Managing the whitelist later  [you]
- **Add or remove someone:** edit the Emails list in this policy. Effective immediately,
  no redeploy.
- **Sharing one list across apps:** if you later host more than one private tool, make the
  list a reusable **Access Group** (Zero Trust → Access → Groups) and reference the group
  from each app's policy, so you maintain it in one place.
- **As code (optional):** the list can be managed via the Cloudflare API or the
  `cloudflare_zero_trust_access_*` Terraform resources. Keep any such file in a **private**
  location, never in this repo — see the warning below.

> **Do not commit the email list into the portal folder.** Everything in the deployed
> directory is served as a static file, so an address list placed there would be readable
> by anyone who passes Access — and, before Part B is finished, by the public. The
> whitelist belongs in the Cloudflare dashboard (or a private Terraform repo), never under
> the served tree. `test/check_integrity.py` fails if it finds a collaborator list in the
> served folder, as a backstop.

### B5. Test the gate  [you]
Open the domain in a **private/incognito** window. You should be sent to Google (or the
PIN screen) before the page loads, and be **refused** if you authenticate with an address
not on the list. Then log in as yourself and confirm the portal loads and works.

--------------------------------------------------------------------------------
## Updating the site later

With the repo connected, **every push redeploys automatically** — no dashboard step. Flow:

1. If the model changed, rebuild assets:
   `python scripts/build_web_assets.py --version vX.Y --sdrd .../SDRD --src ... --ref-root ...`
2. **Run the full test suite, including the NL4-3 round-trip.** Do not deploy on a
   failing round-trip.
3. Commit and push. Cloudflare rebuilds in seconds.
4. Cut a new model as a new `data/vX.Y/` directory; never overwrite a released one, so a
   citation keeps resolving (`docs/VERSIONING.md`).

The immutable-release caching in `_headers` means a browser may hold an old
`data/vX.Y/` copy; that is correct, because a release directory never changes. Code and
`releases.json` are sent `no-cache`, so new versions and UI changes appear on reload.

--------------------------------------------------------------------------------
## What this does and does not protect

- **Protects:** who can open the page. The unpublished atlas and the tool are not on the
  public web; only your listed collaborators get in.
- **Does not need to protect the sequences**, because they never leave the browser — that
  is the static architecture, unchanged by hosting. Access authenticating a user does not
  give it the sequences; there is nothing to give.
- **Logs:** Cloudflare records authentication events (who logged in, when) — access
  metadata, not sequence data, ~24 h retention on free. If even that is sensitive, host
  behind your institution's SSO instead (`docs/HOSTING.md`).
- If collaborators will paste **real patient-derived** sequences, confirm the hosting
  choice with your IRB / IBC / data-governance office even though no data is transmitted.

--------------------------------------------------------------------------------
## Optional: deploy from CI instead of the dashboard Git integration

The dashboard Git integration (Part A) is simplest and needs no secrets. If you prefer
CI, `.github/workflows/deploy-cloudflare.yml` is included; it needs two repository secrets
(`CLOUDFLARE_API_TOKEN` with the Pages:Edit permission, and `CLOUDFLARE_ACCOUNT_ID`) and
uses `wrangler.toml`. Use one path or the other, not both.
