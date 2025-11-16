data(ekman, package = "smacof")
ekman <- 1 - ekman
ekmanData <- makeMDSData(ekman, ekman ^ 2)
ekmanLabels <- as.character(attr(ekman, "Labels"))

