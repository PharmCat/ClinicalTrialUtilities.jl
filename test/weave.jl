using Weave, ClinicalTrialUtilities
weave(joinpath(dirname(pathof(ClinicalTrialUtilities)), "..", "test", "validation_report.jmd");
doctype = "pandoc2pdf",
out_path = :pwd,
pandoc_options=["--toc", "-V colorlinks=true" , "-V linkcolor=blue", "-V urlcolor=red", "-V toccolor=gray"])
