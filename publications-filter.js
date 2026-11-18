// CNA-handling filter for the "List of publications" table in
// subclonal_deconvolution_guide.md. Only that table is affected.
document.addEventListener("DOMContentLoaded", function () {
  // Pandoc does not carry the caption's #id onto the <table>, so find it by caption.
  var table = Array.from(document.querySelectorAll("table")).find(function (t) {
    return t.caption && t.caption.textContent.trim() === "List of publications";
  });
  var select = document.getElementById("publications-cna-filter");
  var count = document.getElementById("publications-count");
  if (!table || !select) return;

  var rows = Array.from(table.tBodies[0].rows);

  // CNA handling is column 5 (index 4). Use the text before any " (" qualifier,
  // so "Yes (via ω)" counts as "Yes".
  var category = function (row) {
    return row.cells[4].textContent.split(" (")[0].trim();
  };

  var apply = function () {
    var value = select.value, shown = 0;
    rows.forEach(function (row) {
      var show = value === "" || category(row) === value;
      row.style.display = show ? "" : "none";
      if (show) shown++;
    });
    count.textContent = "Showing " + shown + " of " + rows.length;
  };

  select.addEventListener("change", apply);
  apply();
});
