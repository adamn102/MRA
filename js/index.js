// SPDX-License-Identifier: MIT
import { bootstrapPage, privacyPanel, note, NEF_NOTE } from "./page.js";

const page = await bootstrapPage("index.html");
document.getElementById("privacy-slot").appendChild(privacyPanel());

const region = page.release.modelled_region;
document.getElementById("scope-slot").appendChild(note(
  `The portal covers <strong>NL4-3 reference positions ${region.first_position.toLocaleString()}` +
  `&ndash;${region.last_position.toLocaleString()}</strong>, where position 1 is the gag start ` +
  `codon. That is the extent of the fitted model. ` + NEF_NOTE
));

// Every link carries its version explicitly, so a copied URL keeps resolving to this release.
const v = page.version.replace(/^v/, "");
for (const [id, href] of [["link-browse", "browse.html"], ["link-sequence", "sequence.html"],
                          ["link-batch", "batch.html"], ["link-download", "download.html"]]) {
  document.getElementById(id).href = href + "?v=" + v;
}
