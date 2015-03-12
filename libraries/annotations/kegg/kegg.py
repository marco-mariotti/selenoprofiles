#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/04/19 13:41:10
#
# This script 

import sys,os,re
from SOAPpy import *
from urllib import urlretrieve


def main():
    pw = 'hsa:04810'
    genesRed  = ['hsa:5594','hsa:2149','hsa:4659','hsa:10672','hsa:1398']
    genesBlue = ['hsa:2335','hsa:8503','hsa:4638','hsa:2909','hsa:5605','hsa:93408','hsa:5747','hsa:673','hsa:624','hsa:91807','hsa:1131','hsa:26230','hsa:4342','hsa:5880','hsa:5500','hsa:5595','hsa:1399','hsa:3655','hsa:1129','hsa:10398','hsa:5216','hsa:1133','hsa:6654','hsa:9564','hsa:5291','hsa:23533','hsa:369','hsa:5604','hsa:2147','hsa:5293','hsa:3685','hsa:3694','hsa:3687','hsa:5217','hsa:5894','hsa:6655','hsa:56924','hsa:3695','hsa:375189','hsa:89846','hsa:7410','hsa:1128','hsa:23365','hsa:2768','hsa:5296','hsa:5290','hsa:1445','hsa:3675','hsa:998','hsa:723961','hsa:3985','hsa:4636','hsa:3693','hsa:5295','hsa:22800','hsa:6237','hsa:22808','hsa:57144','hsa:7409','hsa:10451','hsa:4893','hsa:3680','hsa:3673','hsa:5294','hsa:387','hsa:5063','hsa:5501','hsa:3845','hsa:4633','hsa:103910','hsa:3674','hsa:3691','hsa:345456','hsa:5879','hsa:10627','hsa:58498','hsa:29895','hsa:623','hsa:6093','hsa:1072']
    genes     = ['hsa:2147','hsa:2149','hsa:929','hsa:10672','hsa:2768','hsa:9138','hsa:23365','hsa:387','hsa:6093','hsa:9475','hsa:1729','hsa:1730','hsa:5216','hsa:5217','hsa:345456','hsa:375189','hsa:4659','hsa:5499','hsa:5500','hsa:5501','hsa:4638','hsa:85366','hsa:91807','hsa:103910','hsa:10398','hsa:10627','hsa:29895','hsa:4633','hsa:4636','hsa:58498','hsa:93408','hsa:1950','hsa:1956','hsa:3645','hsa:723961','hsa:6654','hsa:6655','hsa:6237','hsa:3845','hsa:22800','hsa:22808','hsa:4893','hsa:5291','hsa:5294','hsa:5290','hsa:23533','hsa:5293','hsa:5295','hsa:5296','hsa:8503','hsa:10451','hsa:7410','hsa:7409','hsa:26230','hsa:7074','hsa:369','hsa:673','hsa:4342','hsa:5894','hsa:5604','hsa:5605','hsa:5594','hsa:5595','hsa:5880','hsa:5879','hsa:8874','hsa:9459','hsa:5063','hsa:57144','hsa:10298','hsa:5058','hsa:5062','hsa:56924','hsa:3984','hsa:3985','hsa:1072','hsa:1073','hsa:2909','hsa:2335','hsa:3683','hsa:3675','hsa:3674','hsa:8516','hsa:3693','hsa:3655','hsa:3696','hsa:3694','hsa:3691','hsa:3679','hsa:22801','hsa:3685','hsa:3695','hsa:3687','hsa:8515','hsa:3680','hsa:3688','hsa:3681','hsa:3678','hsa:3684','hsa:3672','hsa:3673','hsa:5747','hsa:9564','hsa:1399','hsa:1398','hsa:1445','hsa:1128','hsa:1129','hsa:1131','hsa:1133','hsa:623','hsa:624','hsa:55970','hsa:2245','hsa:89846','hsa:998']
    PSGs      = ['hsa:2260', 'hsa:9138', 'hsa:6655', 'hsa:7410', 'hsa:1729', 'hsa:2245', 'hsa:1956', 'hsa:10298', 'hsa:3684']


class Kegg:
    """
    KEGG object take in argument pw id e.g.: 'hsa:04810'.
    takes also list(s) of genes (KEGG IDs), colors to paint them, and if we want
    to paint background or forground.
    
    Examples:
    *********
    # retrieve image of one pathway
    k = Kegg('hsa:04810')
    k.getImg('lala.png')

    # retrive image of one pathway coloring in green a list of genes, open in firefox:
    k = Kegg('hsa:04810',Ggenes = ['hsa:5594', 'hsa:2149', 'hsa:4659', 'hsa:10672', 'hsa:1398'],colors = 'green')
    os.system('firefox '+k.pngurl)

    # retrieve image of one pathway coloring in blue a list of genes, in red the other,
    in purple the junction (last color of the list), and in yellow letter a third list, open in firefox:
    k = Kegg('hsa:04810',Ggenes = [['hsa:5594', 'hsa:2149', 'hsa:4659', 'hsa:10672', 'hsa:1398'],\
    ['hsa:5879', 'hsa:10627', 'hsa:58498', 'hsa:29895', 'hsa:623', 'hsa:6093', 'hsa:1072'],\
    ['hsa:2260', 'hsa:9138', 'hsa:6655', 'hsa:7410', 'hsa:1729', 'hsa:2245', 'hsa:1956', 'hsa:10298', 'hsa:3684']],\
    colors = ['red','blue','yellow','purple'],bg = [1,1,0])
    os.system('firefox '+k.pngurl)
    """

    def __init__(self, pw, Ggenes = [[]], colors = \
                 ['grey','pink', 'cyan','red','orange', \
                  'yellow', 'green', 'blue', 'indigo', 'violet'],\
                 bg=[1]):
        '''
        You can either path lists of colors, either just one color
        '''
        
        if list(Ggenes[0])!=Ggenes[0]:
            Ggenes = [Ggenes]
        try: list(bg)
        except TypeError: bg = [bg]
        if list(colors) != colors:
            colors = [colors]
        wsdl = 'http://soap.genome.jp/KEGG.wsdl'
        if len (Ggenes) > len (colors):
            sys.exit("ERROR: not enough colors given")
        if bg == []:
            bg = [True]*len(Ggenes)
        self.pw = pw
        self.Ggenes   = Ggenes
        self.bg       = bg
        self.serv     = WSDL.Proxy(wsdl)
        self.elements = self.serv.get_elements_by_pathway(\
            re.sub('([a-z]{3}):', r'path:\1', self.pw))
        self.bgcolors = ['white']*len(self)
        self.fgcolors = ['black']*len(self)
        self.pngurl   = self._paintPW(colors)

    def __len__(self):
        '''
        returns number of elements in pathway
        '''
        return len(self.elements)

    def _paintPW(self, colors):
        '''
        return String with url of png
        '''
        Elts = map (lambda x: x.element_id, self.elements)
        print Elts
        for e in Elts:
            e = e - 1
            if not self.elements[e].type == 'gene': continue
            for G in range(0, len(self.Ggenes)):
                if len (set(self.elements[e].names) \
                        & set(self.Ggenes[G]))==0: continue
                if self.bgcolors[e] != 'white' \
                       and self.bgcolors[e] != 'grey' \
                       and self.bgcolors[e] != colors[G] \
                       and self.bg[G] != 0:
                    if self.bg[G]: self.bgcolors[e] = colors[-1]
                    else: self.fgcolors[e] = colors[G]
                else:
                    if self.bg[G]: self.bgcolors[e] = colors[G]
                    else: self.fgcolors[e] = colors[G]
        return self.serv.color_pathway_by_elements(\
            re.sub('hsa:','path:hsa',self.pw), Elts,\
            self.fgcolors,\
            self.bgcolors)
    def getImg(self, outfile):
        '''
        retrieve png from url, to file.
        '''
        urlretrieve(self.pngurl, outfile)


if __name__ == "__main__":
    sys.exit(main())

