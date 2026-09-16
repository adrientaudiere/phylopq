#' Objects re-exported from MiscMetabar
#'
#' @description
#' `align_pq()` and `is_mafft_installed()` were introduced in phylopq 0.2.0 and
#' now live in MiscMetabar, so that every package of the pqverse can reach the
#' MAFFT backend without depending on phylopq — MiscMetabar must keep working
#' on its own, so it cannot depend on phylopq, and it is the one package the
#' others already share.
#'
#' They are re-exported here, so `phylopq::align_pq()` and
#' `phylopq::is_mafft_installed()` keep working unchanged. The only visible
#' difference is the option that holds the path to the executable, renamed from
#' `phylopq.mafftpath` to `MiscMetabar.mafftpath`.
#'
#' @name reexports
#' @keywords internal
#' @seealso [MiscMetabar::align_pq()], [MiscMetabar::is_mafft_installed()],
#'   [delim_pq()], [delim_multi_pq()]
NULL

#' @importFrom MiscMetabar align_pq
#' @export
MiscMetabar::align_pq

#' @importFrom MiscMetabar is_mafft_installed
#' @export
MiscMetabar::is_mafft_installed
