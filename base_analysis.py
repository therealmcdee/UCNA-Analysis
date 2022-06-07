import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas
import sys
import scipy.optimize as spyopt 
import scipy.special
#----------------------------------------------------------------------------
def Weights(x,locs):
    wghts = []
    alen = 360. 

    for i in range(0,256,2):
        vx0 = locs[i] - x[0]
        vy0 = locs[i+1] -x[1]
        if i < 254:
            vx1 = locs[i+2]-x[0]
            vy1 = locs[i+3]-x[1]
        else:
            vx1 = locs[0] - x[0]
            vy1 = locs[1] - x[1]
        mag1 = np.sqrt(vx0*vx0 + vy0*vy0)
        mag2 = np.sqrt(vx1*vx1 + vy1*vy1)

        v1dotv2 = (vx0*vx1 + vy0*vy1)/(mag1*mag2)
        aveg = (mag1+mag2)/2
        v1dotv2 = x[2]*np.arccos(v1dotv2)/(2*np.pi)*np.exp(-aveg/alen)
        wghts.append(v1dotv2)

    wghts = np.asarray(wghts)
    hits = np.zeros(128)
    df = 0

    for i in range(0,128):
        df = df + np.power((wghts[i] - locs[256+i]),2)

    df = df/128

    return df


#----------------------------------------------------------------------------
def build_map():
    x = []
    y = []

    read = open('data_files/sipm_map.txt','r')
    for line in read:
        seg = line.split()
        x.append(float(seg[1]))
        y.append(float(seg[2]))
    read.close()
    x = np.asarray(x)
    y = np.asarray(y)

    return x,y

def calc_sipm_pos():
    pos = np.zeros(128*3)
    dx = 0.429
    theta = 11.25
    rad = 8.794
    for i in range(0,16):
        comp_ang = (theta*(1+2*i)/180+.5)*np.pi
        for j in range(0,8):
            x = rad*np.cos(i*2*theta*np.pi/180)+dx*j*np.cos(comp_ang)
            y = rad*np.sin(i*2*theta*np.pi/180)+dx*j*np.sin(comp_ang)
            pos[(i*8+j)*2] = x
            pos[(i*8+j)*2+1] = y
    return pos

#----------------------------------------------------------------------------

def main(argv):

    nstart = 8

    sipm_x, sipm_y = build_map()
    sipm_locs = calc_sipm_pos()
    sipm_id = []
    #three boards; 24 'working' channels
    for h in range(3,0,-1):
        for j in range(0,8):
            sipm_id.append(j%8+h*8)

    nrun = 12
    if len(argv) > 1 :
        nrun = int(argv[1])

    # open the data as a pandas data frame
    df = pandas.read_csv('data_files/{}.data'.format(nrun),delimiter=";")
    
    # create a 16 fold histogram array to hold the individual channel data.
    LG_Histos = [[0 for x in range(0,1024)] for y in range(0,24)]
    TT_Histos = []
    evt_HG = []
    time = []
    weighted_pos_x = []
    weighted_pos_y = []
    resx = []
    resy = []
    events = 0

    for index,row in df.iterrows():
        for nE in range(0,int(row['NEventsInCluster'])) :
            nch = int(row['ASIC_{}'.format(nE)])
            if nch== 0: 
                total_charge = 0
                adc_array = []
                num_x = 0
                num_y = 0
                for j in range(nstart,nstart+24):
                    adc_value = int(row['CHARGE_HG_{}_{}'.format(nch,j)]/16)
                    adc_array.append(adc_value) 
                    if adc_value > 0 :
                        LG_Histos[j-nstart][int(adc_value)] = LG_Histos[j-nstart][int(adc_value)] + 1
                        total_charge = total_charge + adc_value
                if total_charge > 0:
                    evt_HG.append(total_charge)
                for p in range(0,len(sipm_id)):
                    num_x += sipm_x[p]*adc_array[sipm_id[p]-8]
                    num_y += sipm_y[p]*adc_array[sipm_id[p]-8]
                    sipm_locs[p+256] = adc_array[sipm_id[p]-8]
                weighted_pos_x.append(np.sum(num_x)/np.sum(adc_array))
                weighted_pos_y.append(np.sum(num_y)/np.sum(adc_array))

                res = spyopt.minimize(Weights,x0=(2*weighted_pos_x[events],2*weighted_pos_y[events],evt_HG[events]),args=(sipm_locs),method='CG')
                resx.append(res.x[0])
                resy.append(res.x[1]) 
                events += 1
        time.append(float(row['RUN_EventTimecode_ns_0'])*1e-9*0.5/25)

    print('Runtime = {} seconds'.format(time[-1]))
    print('Total Events = {}'.format(events))

    # make a x-axis for the channel histograms
    xhist = np.arange(0,1024,1)
    # create a 4x4 array of histogram for the channels
    fig,axs = plt.subplots(6,4)

    for i in np.arange(0,6):
        for j in np.arange(0,4):
            axs[i][j].step(xhist,LG_Histos[i*4 + j],label='Chn {}'.format(i*4+j))
            axs[i][j].set_xlabel("ADC [n]")
            axs[i][j].set_ylabel("Counts")

    # plot the summed HG ADC
    fig2 = plt.figure()
    plt.hist(evt_HG,bins=1024)

    fig3 = plt.figure()
    plt.hist(time,bins=int(time[-1])+20)#,range=(0,int(local_time[-1]+20)))
    plt.xlabel('Seconds')
    plt.ylabel('Number of Events')

    fig4 = plt.figure()
    plt.scatter(weighted_pos_x,weighted_pos_y,color='b')
    plt.scatter(resx,resy,color='r')
    poly_x = np.zeros(17)
    poly_y=np.zeros(17)
    for i in range(17):
        poly_x[i] = 8.626*np.cos(2*np.pi*i/16)
        poly_y[i] = 8.626*np.sin(2*np.pi*i/16)
    plt.plot(poly_x,poly_y,color='k')
    plt.xlabel(r'$x_{recon}$')
    plt.ylabel(r'$y_{recon}$')

    plt.show()



if __name__=="__main__":
    main(sys.argv)
