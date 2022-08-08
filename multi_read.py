import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas
import sys


def main(argv):

    nrun = []
    evt_HG = []
    if len(argv) > 1 :
        print(argv)
        for i in range(1,len(argv)):
            nrun.append(int(argv[i]))
            evt_HG.append([])

    for f in range(len(nrun)):

    # open the data as a pandas data frame
        df = pandas.read_csv('data_files/{}.data'.format(nrun[f]),delimiter=";")
    
    # create a 16 fold histogram array to hold the individual channel data.
        LG_Histos = [[0 for x in range(0,1024)] for y in range(0,24)]
        TT_Histos = []

        nstart = 8

        for index,row in df.iterrows():
            for nE in range(0,int(row['NEventsInCluster'])) :
                nch = int(row['ASIC_{}'.format(nE)])
                if nch== 0: 
                    total_charge = 0
                    for j in range(nstart,nstart+24):
                        adc_value = int(row['CHARGE_HG_{}_{}'.format(nch,j)]/16)
                        if adc_value > 0 :
                            LG_Histos[j-nstart][int(adc_value)] = LG_Histos[j-nstart][int(adc_value)] + 1
                            total_charge = total_charge + adc_value
                    if total_charge > 0:
                        evt_HG[f].append(total_charge)
                        

    # make a x-axis for the channel histograms
    #xhist = np.arange(0,1024,1)
    # create a 4x4 array of histogram for the channels
    #fig,axs = plt.subplots(6,4)

    #for i in np.arange(0,6):
        #for j in np.arange(0,4):
            #axs[i][j].step(xhist,LG_Histos[i*4 + j],label='Chn {}'.format(i*4+j))
            #axs[i][j].set_xlabel("ADC [n]")
            #axs[i][j].set_ylabel("Counts")


    # plot the summed HG ADC
    posers = ['Pos 1','Pos 2','Pos 3','Pos 4','Pos 5']

    fig2 = plt.figure()
    for i in range(len(evt_HG)):
        plt.hist(evt_HG[i],bins=1024,histtype='step',label='{}'.format(posers[i]))#, ADC Overall = {}'.format(posers[i],np.sum(evt_HG[i])))
    plt.xlabel('ADC',fontsize=18)
    plt.ylabel('Number of Events',fontsize=18)
    plt.title('400 nm, 1 kHz LED',fontsize = 20) 
    plt.tick_params(labelsize=14)
    plt.legend()
    plt.show()



if __name__=="__main__":
    main(sys.argv)
