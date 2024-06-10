lines <- c("time_units: generations",
           "defaults:",
           "  epoch: {start_size: 500, end_time: 0}",
           "demes:",
           "- name: ancestral",
           "  epochs:",
           "  - {start_size: 10000, end_time: 50}",
           paste("- {name: line",
                 1:100,
                 ", ancestors: [ancestral]}", 
                 sep = ""))

fileConn<-file("demo.yaml")
writeLines(lines, fileConn)
close(fileConn)

