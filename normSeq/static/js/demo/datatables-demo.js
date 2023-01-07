// Call the dataTables jQuery plugin
$(document).ready(function() {
  // Query - Table methods
  $('#dataTableMethods').DataTable({
    "order":[[1, "desc"]]
  });
  $('#dataTableEdgeR').DataTable({
    "order":[[4, "asc"]]
  });
  $('#dataTableDeseq').DataTable({
    "order":[[4, "asc"]]
  });
  $('#dataTableNoiseq').DataTable({
    "order":[[3, "asc"]]
  });
  $('#dataTableTtest').DataTable({
    "order":[[3, "asc"]]
  });
  
}); 