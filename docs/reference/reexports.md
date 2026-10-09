# Objects re-exported from MiscMetabar

[`align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html)
and
[`is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mafft_installed.html)
were introduced in phylopq 0.2.0 and now live in MiscMetabar, so that
every package of the pqverse can reach the MAFFT backend without
depending on phylopq — MiscMetabar must keep working on its own, so it
cannot depend on phylopq, and it is the one package the others already
share.

They are re-exported here, so
[`phylopq::align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html)
and
[`phylopq::is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mafft_installed.html)
keep working unchanged. The only visible difference is the option that
holds the path to the executable, renamed from `phylopq.mafftpath` to
`MiscMetabar.mafftpath`.

## See also

[`MiscMetabar::align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html),
[`MiscMetabar::is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mafft_installed.html),
[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md),
[`delim_multi_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_multi_pq.md)
