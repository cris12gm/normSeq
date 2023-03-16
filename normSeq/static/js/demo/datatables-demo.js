// Call the dataTables jQuery plugin
$(document).ready(function() {
  // Query - Table methods
  $('#dataTableMethods').DataTable({
    "order":[[1, "desc"]]
  });
  $('.dataTableEdgeR').DataTable({
    "order":[[5, "asc"]]
  });
  $('.dataTableDeseq').DataTable({
    "order":[[5, "asc"]]
  });
  $('.dataTableNoiseq').DataTable({
    "order":[[4, "asc"]]
  });
  $('.dataTableTtest').DataTable({
    "order":[[4, "asc"]]
  });
  $('.dataTableSummaryDE').DataTable({
    
  });
  $('.dataTableConsensus').DataTable({
    "order":[[1, "asc"]]
  });
  
  
}); 