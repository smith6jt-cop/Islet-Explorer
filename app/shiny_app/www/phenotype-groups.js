// Grouped phenotype choices (R/00_globals.R phenotype_group_checkboxes(); owner, 2026-10-27): the All | None links
// of a group, or of the whole list, check or clear its boxes, then tell Shiny once.
(function () {
function phenotypeBoxes(scope) {
  var name = scope.closest('.pheno-group-checkboxes').attr('id');
  return scope.find('input:checkbox').filter(function () { return this.name === name; });
}

function updateMacrophageToggles(scope) {
  scope.find('.pheno-group-toggle').each(function () {
    var boxes = phenotypeBoxes($(this).closest('.pheno-group'));
    var selected = boxes.filter(':checked').length;
    this.checked = boxes.length > 0 && selected === boxes.length;
    this.indeterminate = selected > 0 && selected < boxes.length;
  });
}

$(document).on('click', '.pheno-select', function (event) {
  event.preventDefault();
  var link = $(this);
  var scope = link.hasClass('pheno-scope-group') ? link.closest('.pheno-group') : link.closest('.pheno-group-checkboxes');
  // only the list's own boxes: a Macrophage heading also holds the marker checkboxes, a separate input
  var boxes = phenotypeBoxes(scope);
  boxes.prop('checked', link.hasClass('pheno-select-all'));
  boxes.first().trigger('change');
  updateMacrophageToggles(scope);
});

// A Macrophage checkbox selects its displayed children without adding a parent
// value to Shiny's phenotype input or changing the four marker selections.
// Click runs before the checkbox's change event reaches Shiny's ancestor
// binding. Update children first, so synchronization cannot restore an old
// indeterminate/unchecked state before this handler reads the new choice.
$(document).on('click', '.pheno-group-toggle', function () {
  var group = $(this).closest('.pheno-group');
  var boxes = phenotypeBoxes(group);
  boxes.prop('checked', this.checked);
  boxes.first().trigger('change');
  updateMacrophageToggles(group);
});

$(document).on('change', '.pheno-group-checkboxes input:checkbox', function () {
  updateMacrophageToggles($(this).closest('.pheno-group-checkboxes'));
});
$(document).on('shiny:bound', function (event) {
  updateMacrophageToggles($(event.target).closest('.pheno-group-checkboxes'));
});
$(function () { updateMacrophageToggles($(document)); });
})();
