from ROOT import TRandom3, TMath
from optparse import OptionParser
import progressbar

__author__ = 'Diego Alejandro Sanz Becerra'

h = 0
t = 1

class bet:
    def __init__(self, N=20):
        self.N = N
        self.rand = TRandom3(123124)
        self.totalHeadsNum = 0
        self.totalTailsNum = 0
        self.headsNplus1 = 0
        self.headsNtail1 = 0
        self.results = {}

    def toss(self):
        #  start counters
        heads = 0
        tails = 0
        #  toss N times
        for i in xrange(self.N):
            toss = TMath.Nint(self.rand.Rndm())
            if toss is h:
                self.totalHeadsNum += 1
                heads += 1
            else:
                self.totalTailsNum += 1
                tails += 1
        #  toss N + 1 event
        toss = TMath.Nint(self.rand.Rndm())
        if toss is h:
            self.totalHeadsNum += 1
            #  if the previous were all heads (this is head too)
            if heads is self.N:
                self.headsNplus1 += 1
        else:
            self.totalTailsNum += 1
            #  if it the previous one were all heads (this is tail)
            if heads is self.N:
                self.headsNtail1 += 1

    def playGame(self, trials = 1000):
        self.trials = trials
        print "STARTING GAME WITH ", trials, " TRIALS:"
        widgets = [
                progressbar.Percentage(),
                ' ', progressbar.Bar(marker='>'),
                ' ', progressbar.Timer(),
                ' ', progressbar.ETA()  #, DA: this two work great!
                # ' ', progressbar.AdaptiveETA(),
                # ' ', progressbar.AdaptiveTransferSpeed(),
                ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=trials)
        bar.start()
        for i in xrange(trials):
            self.toss()
            bar.update(i + 1)
        bar.finish()
        print "RESULTS:"

    def probabilityHeadsNtails1(self):
        return float(self.headsNtail1)/float(self.headsNtail1 + self.headsNplus1)

    def probabilityHeadsNplus1(self):
        return float(self.headsNplus1)/float(self.headsNtail1 + self.headsNplus1)

    def probabilityHeads(self):
        return float(self.totalHeadsNum)/float(self.totalHeadsNum + self.totalTailsNum)

    def probabilityTails(self):
        return float(self.totalTailsNum)/float(self.totalHeadsNum + self.totalTailsNum)

    def get_results(self):
        if (self.headsNplus1 + self.headsNtail1) is 0:
            print "Not enough trials."
        else:
            self.results['Nheads1tail'] = self.probabilityHeadsNtails1()
            self.results['Nplus1heads'] = self.probabilityHeadsNplus1()
            self.results['heads'] = self.probabilityHeads()
            self.results['tails'] = self.probabilityTails()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-n', '--en', dest='n', default=5, type='int', help='Number of consecutive heads: e.g. 20')
    parser.add_option('-t', '--trials', dest='t', default=10000000, type='int', help='Number of trials to test the bet: e.g. 1000')
    (options, args) = parser.parse_args()
    n = int(options.n)
    trials = int(options.t)
    z = bet(n)
    z.playGame(trials)
    z.get_results()
    print z.results
