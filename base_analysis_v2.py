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

def Gaussian(x,a,b,c):
    return a*np.exp(-np.power((x-b),2)/(c*c*2))

def calc_avg(data):
    return np.sum(data)/len(data)

def calc_stddev(data):
    data_avg = np.sum(data)/float(len(data))
    num = []
    for i in range(len(data)):
        num.append((data[i]-data_avg)**2)
    sig = np.sqrt(np.sum(num)/(len(data)-1))
    return sig

def find_list_max(mylist):
    maxi = mylist[0]
    index = 0
    for i in range(len(mylist)):
        if mylist[i] > maxi:
            maxi = mylist[i]
            index = i
    return index

#----------------------------------------------------------------------------

def main(argv):

    nstart = 0

    sipm_x, sipm_y = build_map()
    sipm_locs = calc_sipm_pos()
    sipm_id = []
    #three boards; 24 'working' channels
    for p in range(3):
        for h in range(31,-1,-1):
            sipm_id.append(h)

    nrun = 12
    if len(argv) > 1 :
        nrun = int(argv[1])

    # open the data as a pandas data frame
    df = pandas.read_csv('data_files/{}.data'.format(nrun),nrows = 100,delimiter=";")
    
    # create a 16 fold histogram array to hold the individual channel data.
    LG_Histos = [[0 for x in range(0,1024)] for y in range(0,96)]
    TT_Histos = []
    evt_HG = []
    single_HG = []
    double_HG = []
    triple_HG = []
    time = []
    time_rollover = 0
    delta_t = []
    time_holder = []
    start_time = 0
    pos_x = []
    pos_y = []
    resx = []
    resy = []
    events = 0
    triple_events = 0
    double_events = 0
    single_events = 0
    ADC_thresh = 60
    mult_array = [0 for x in range(71)]
    mult_array_A = []
    mult_array_B = []
    mult_array_C = []
    single_A = 0
    single_B = 0
    single_C = 0

    for index,row in df.iterrows():
        print('EVENT: {}'.format(index))
        total_charge = 0
        adc_array = np.zeros(96)
        adc_value = 0
        num_x = 0
        num_y = 0
        multiplicity = 0
        mult_A = 0
        mult_B = 0
        mult_C = 0
        num_events = int(row['NEventsInCluster'])
        for n_eve in range(int(row['NEventsInCluster'])):
            nch = int(row['ASIC_{}'.format(n_eve)])
            if nch == 1:
                for j in range(32):
                    adc_array[j] = int(row['CHARGE_HG_{}_{}'.format(n_eve,sipm_id[j])]/16)
                    adc_value = int(row['CHARGE_HG_{}_{}'.format(n_eve,sipm_id[j])]/16)
                    if adc_value > 0:
                        LG_Histos[j][int(adc_value)] = LG_Histos[j][int(adc_value)] + 1
                        total_charge = total_charge + adc_value
                        if adc_value>ADC_thresh:
                            multiplicity += 1
                            mult_B += 1
            elif nch == 0:
                for j in range(32):
                    adc_array[j+32] = int(row['CHARGE_HG_{}_{}'.format(n_eve,sipm_id[j+32])]/16)
                    adc_value = int(row['CHARGE_HG_{}_{}'.format(n_eve,sipm_id[j+32])]/16)
                    if adc_value > 0:
                        LG_Histos[j+32][int(adc_value)] = LG_Histos[j+32][int(adc_value)] + 1
                        total_charge = total_charge + adc_value
                        if adc_value > ADC_thresh:
                            multiplicity += 1
                            mult_A += 1
            elif nch == 2:
                for j in range(32):
                    adc_array[j+64] = int(row['CHARGE_HG_{}_{}'.format(n_eve,sipm_id[j+64])]/16)
                    adc_value = int(row['CHARGE_HG_{}_{}'.format(n_eve,sipm_id[j+64])]/16)
                    if adc_value > 0:
                        LG_Histos[j+64][int(adc_value)] = LG_Histos[j+64][int(adc_value)]+1
                        total_charge = total_charge + adc_value
                        if adc_value>ADC_thresh:
                            multiplicity += 1
                            mult_C += 1
            if num_events == 1 and nch == 1:
                single_B += 1
            elif num_events == 1 and nch == 0:
                single_A += 1
            elif num_events == 1 and nch ==2:
                single_C += 1

        if num_events == 1:
            single_HG.append(total_charge)
        elif num_events == 2:
            double_HG.append(total_charge)
        elif num_events == 3:
            triple_HG.append(total_charge)
        evt_HG.append(total_charge)
        for p in range(96):
            num_x += sipm_x[p]*adc_array[p]
            num_y += sipm_y[p]*adc_array[p]
            sipm_locs[p+256] = adc_array[p]
        pos_x.append(num_x/np.sum(adc_array))
        pos_y.append(num_y/np.sum(adc_array))
        res = spyopt.minimize(Weights,x0=(pos_x[events],pos_y[events],evt_HG[events]),args=(sipm_locs))#,method='CG')
        resx.append(res.x[0])
        resy.append(res.x[1])

        mult_array[multiplicity] = mult_array[multiplicity] + 1
        mult_array_A.append(mult_A)
        mult_array_B.append(mult_B)
        mult_array_C.append(mult_C)

        time_holder.append(float(row['CLUSTER_RUN_Timecode_ns']))
        if events == 0:
            start_time = float(row['CLUSTER_RUN_Timecode_ns'])
        elif events > 0:
            delta_t.append((time_holder[events]-time_holder[events-1])*0.5e-9/25)
        if float(row['CLUSTER_RUN_Timecode_ns']) - start_time < 0:
            if time_rollover == 0:
                time_rollover += 1
                time.append(time[events-1]+float(row['CLUSTER_RUN_Timecode_ns'])*0.5e-9/25)
            else:
                time.append(time[events-1]+float(row['CLUSTER_RUN_Timecode_ns'])*0.5e-9/25)
        else:
            time.append((float(row['CLUSTER_RUN_Timecode_ns'])-start_time)*0.5e-9/25)
        
        events+=1
        if n_eve == 1:
            double_events += 1
        elif n_eve == 2:
            triple_events += 1
        elif n_eve == 0:
            single_events += 1
            
    print('Runtime = {} seconds'.format(time[-1]))
    print('Total Events = {}'.format(double_events+single_events))
    print('Triple Events = {}'.format(triple_events))
    print('Double Events = {}'.format(double_events))
    print('Single Events = {}'.format(single_events))
    print('Single Events Board A = {}'.format(single_A))
    print('Single Events Board B = {}'.format(single_B))

    # make a x-axis for the channel histograms
    xhist = np.arange(0,1024,1)

    ebin = []
    ebin2 = []

    # create a 4x4 array of histogram for the channels
    fig,axs = plt.subplots(4,8,sharex=True,sharey=True)
    plt.suptitle('ASIC 0')
    for i in np.arange(0,4):
        for j in np.arange(0,8):
            axs[i][j].step(xhist,LG_Histos[i*8+j+32],label='Chn {}'.format(i*8+j+32))
            #popt,pcov = spyopt.curve_fit(Gaussian,xhist,LG_Histos[i*8+j+32],p0=(find_list_max(LG_Histos[i*8+j+32]),calc_avg(LG_Histos[i*8+j+32]),calc_stddev(LG_Histos[i*8+j+32])))
            #axs[i][j].plot(xhist,Gaussian(xhist,*popt),'--',color='r',label='ADC = {:4.3f}'.format(popt[1]))
            axs[i][j].set_xlabel("ADC [n]")
            axs[i][j].set_ylabel("Counts")

    fig2,axs = plt.subplots(4,8,sharex=True,sharey=True)
    plt.suptitle('ASIC 1')
    for i in np.arange(0,4):
        for j in np.arange(0,8):
            axs[i][j].step(xhist,LG_Histos[i*8+j],label='Chn {}'.format(i*8+j))
            #popt2,pcov2 = spyopt.curve_fit(Gaussian,xhist,LG_Histos[i*8+j],p0=(mAx,calc_avg(LG_Histos[i*8+j]),calc_stddev(LG_Histos[i*8+j])))
            #axs[i][j].plot(xhist,Gaussian(xhist,*popt2),'--',color='r',label='ADC = {:4.3f}'.format(popt2[1]))
            axs[i][j].set_xlabel('ADC [n]')
            axs[i][j].set_ylabel('Counts')
    plt.legend()

    # plot the summed HG ADC
    fig3 = plt.figure()
    plt.title('Total HG ADC')
    plt.hist(evt_HG,bins=1024)

    fig4 = plt.figure()
    plt.title('Total Charge (Single v Double)')
    npe,bpe,pp = plt.hist(single_HG,bins=1024,color='r',label='Single Trigger')
    for i in range(1,len(bpe)):
        ebin.append((bpe[i]+bpe[i-1])/2)
    poptpe,pcovpe = spyopt.curve_fit(Gaussian,ebin,npe,p0=(npe.max(),calc_avg(single_HG),calc_stddev(single_HG)))
    plt.plot(bpe,Gaussian(bpe,*poptpe),'--',color='darkred',label=r'$\langle ADC \rangle$ = {:4.3f}, $\sigma =$ {:4.3f}'.format(poptpe[1],abs(poptpe[2])))

    npe2,bpe2,pp2 = plt.hist(double_HG,bins=1024,color='b',label='Double Trigger')
    for i in range(1,len(bpe2)):
        ebin2.append((bpe2[i]+bpe2[i-1])/2)
    poptpe2,pcovpe2 = spyopt.curve_fit(Gaussian,ebin2,npe2,p0=(npe2.max(),calc_avg(double_HG),calc_stddev(double_HG)))
    plt.plot(bpe2,Gaussian(bpe2,*poptpe2),'--',color='teal',label=r'$\langle ADC \rangle$ = {:4.3f}, $\sigma =$ {:4.3f}'.format(poptpe2[1],abs(poptpe2[2])))
    plt.legend()

    fig5,axs = plt.subplots()
    axs1 = plt.subplot(121)
    axs1.hist(time,bins=1)#int(time[-1]))#,range=(0,int(local_time[-1]+20)))
    axs1.set_xlabel('Seconds')
    axs1.set_ylabel('Number of Events')
    axs2 = plt.subplot(122)
    axs2.hist(delta_t,bins=100)
    axs2.set_xlabel(r'$\delta$t [s]')
    axs2.set_ylabel('Counts')

    fig6 = plt.figure()
    plt.scatter(pos_x,pos_y,color='b',label='ADC Weighted Position')
    plt.scatter(resx,resy,color='lime',label='Fitted Position')
    poly_x = np.zeros(17)
    poly_y=np.zeros(17)
    for i in range(17):
        poly_x[i] = 8.626*np.cos(2*np.pi*i/16)
        poly_y[i] = 8.626*np.sin(2*np.pi*i/16)
    plt.plot(poly_x,poly_y,color='k')
    plt.title('Position Reconstruction (64 Active Channels)',fontsize=20)
    plt.xlabel(r'$x_{recon}$',fontsize=18)
    plt.ylabel(r'$y_{recon}$',fontsize=18)
    plt.legend(loc='upper right')

    mult_x = np.arange(71)

    fig7 = plt.figure()#,axis = plt.subplots(1,2)
    plt.step(mult_x,mult_array)
    plt.title('Total Multiplicity of SiPMs (ADC > {})'.format(ADC_thresh),fontsize=20)
    plt.xlabel('# of SiPMs triggered',fontsize =18)
    plt.ylabel('Counts',fontsize=18)

    #plt.hist2d(mult_array_A,mult_array_B,bins=50)
    #plt.title('SiPMs Triggered: ASIC A v. ASIC B (ADC > {})'.format(ADC_thresh),fontsize=18)
    #plt.xlabel('ASIC A SiPMs Triggered',fontsize=18)
    #plt.ylabel('ASIC B SiPMs Triggered',fontsize=18)

    plt.show()



if __name__=="__main__":
    main(sys.argv)
