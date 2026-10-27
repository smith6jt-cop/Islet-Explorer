// Grouped phenotype choices (R/00_globals.R phenotype_group_checkboxes(); owner, 2026-10-27): the All | None links
// of a group, or of the whole list, check or clear its boxes, then tell Shiny once.
$(document).on('click', '.pheno-select', function (event) {
  event.preventDefault();
  var link = $(this);
  var scope = link.hasClass('pheno-scope-group') ? link.closest('.pheno-group') : link.closest('.pheno-group-checkboxes');
  var boxes = scope.find('input:checkbox');
  boxes.prop('checked', link.hasClass('pheno-select-all'));
  boxes.first().trigger('change');
});
