data(morse, package = "smacof")
morseData <- makeMDSData(morse, 1 / morse)
morseLabels <- as.character(attr(morse, "Labels"))

