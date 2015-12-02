from AbstractClasses.RunClass import Run
from AbstractClasses.newAnalysis import Analysis
from AbstractClasses.AnalysisCollection import AnalysisCollection


runnr = 211

run = Run(211)

ana = Analysis(run)

col = AnalysisCollection()
col.AddAnalysis(ana)

col.ShowSignalVSRate()

raw_input("done.")