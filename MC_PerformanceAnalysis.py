from AbstractClasses.MCPerformance import MCPerformance

MCAnalysis = MCPerformance(verbose=True)
MCAnalysis.DoSignalHeightScan(hits_per_height=300000)

