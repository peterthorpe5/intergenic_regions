"use strict";
document.querySelectorAll("input[data-table]").forEach((input) => {
  input.addEventListener("input", () => {
    const table = document.getElementById(input.dataset.table);
    const query = input.value.toLocaleLowerCase("en-GB").trim();
    let visible = 0;
    table.querySelectorAll("tbody tr").forEach((row) => {
      row.hidden = !row.textContent.toLocaleLowerCase("en-GB").includes(query);
      if (!row.hidden) visible += 1;
    });
    table.closest("section").querySelector(".table-status").textContent =
      `${visible} matching preview rows. Full results are in the TSV file.`;
  });
});
document.querySelectorAll("button[data-sort]").forEach((button) => {
  button.addEventListener("click", () => {
    const column = Number(button.dataset.sort);
    const table = button.closest("table");
    const body = table.querySelector("tbody");
    const direction = button.dataset.order === "asc" ? -1 : 1;
    button.dataset.order = direction === 1 ? "asc" : "desc";
    const rows = Array.from(body.rows);
    rows.sort((a, b) => {
      const left = a.cells[column].textContent.trim();
      const right = b.cells[column].textContent.trim();
      if (left === "NA" || right === "NA") {
        return left === right ? 0 : left === "NA" ? 1 : -1;
      }
      const numeric = left !== "" && right !== "" &&
        Number.isFinite(Number(left)) && Number.isFinite(Number(right));
      return direction * (numeric ? Number(left) - Number(right) :
        left.localeCompare(right, "en-GB", {numeric: true}));
    });
    rows.forEach((row) => body.appendChild(row));
    table.querySelectorAll("th").forEach((heading) => heading.removeAttribute("aria-sort"));
    button.closest("th").setAttribute("aria-sort", direction === 1 ? "ascending" : "descending");
  });
});
