import fitsio
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg') # Don't use X-server.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

plt.style.use('/global/u2/m/mjarvis/SVA1StyleSheet.mplstyle')

tag = ''
#tag = '_good_seeing'  # Uncomment to do "good seeing" run
tag = '_psfex'        # Uncomment to do psfex run

def pretty_rho1(meanr, rho, sig, rho3=None, sig3=None, rho4=None, sig4=None):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho1_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho3 is not None:
        plt.plot(meanr*1.03, rho3, color='green')
        plt.plot(meanr*1.03, -rho3, color='green', ls=':')
        plt.errorbar(meanr[rho3>0]*1.03, rho3[rho3>0], yerr=sig3[rho3>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho3<0]*1.03, -rho3[rho3<0], yerr=sig3[rho3<0], color='green', ls='', marker='s')
        rho3_line = plt.errorbar(-meanr, rho3, yerr=sig3, color='green', marker='s')
    if rho4 is not None:
        plt.plot(meanr*1.06, rho4, color='red')
        plt.plot(meanr*1.06, -rho4, color='red', ls=':')
        plt.errorbar(meanr[rho4>0]*1.06, rho4[rho4>0], yerr=sig4[rho4>0], color='red', ls='', marker='^')
        plt.errorbar(meanr[rho4<0]*1.06, -rho4[rho4<0], yerr=sig4[rho4<0], color='red', ls='', marker='^')
        rho4_line = plt.errorbar(-meanr, rho4, yerr=sig4, color='red', marker='^')

    if rho3 is not None and rho4 is not None:
        plt.legend([rho1_line, rho3_line, rho4_line],
                   [r'$\rho_1(\theta)$', r'$\rho_3(\theta)$', r'$\rho_4(\theta)$'],
                   loc='upper right', fontsize=24)
    else:
        plt.legend([rho1_line],
                   [r'$\rho_1(\theta)$'],
                   loc='upper right')
    plt.ylim( [1.e-12, 1.e-5] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig, rho5=None, sig5=None):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho2_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho5 is not None:
        plt.plot(meanr*1.03, rho5, color='green')
        plt.plot(meanr*1.03, -rho5, color='green', ls=':')
        plt.errorbar(meanr[rho5>0]*1.03, rho5[rho5>0], yerr=sig5[rho5>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho5<0]*1.03, -rho5[rho5<0], yerr=sig5[rho5<0], color='green', ls='', marker='s')
        rho5_line = plt.errorbar(-meanr, rho5, yerr=sig5, color='green', marker='s')

    if rho5 is not None:
        plt.legend([rho2_line, rho5_line],
                   [r'$\rho_2(\theta)$', r'$\rho_5(\theta)$'],
                   loc='upper right', fontsize=24)
    else:
        plt.legend([rho2_line],
                   [r'$\rho_2(\theta)$'],
                   loc='upper right')
    plt.ylim( [1.e-12, 1.e-5] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho0(meanr, rho, sig):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')

    plt.legend([rho0_line],
               [r'$\rho_0(\theta)$'],
               loc='upper right')
    plt.ylim( [1.e-12, 1.e-5] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

file_name = 'run_rho' + tag + '.out'
print('reading ',file_name)
with open(file_name, 'rb') as f:
    results = pickle.load(f)

print('results keys = ',results.keys())

for band in ['u','g','r','i','z','y','riz','all']:
    if (band,'ecat','ecat') not in results:
        continue

    print('Work on band ',band)
    rho0 = results[(band,'ecat','ecat')]
    rho1 = results[(band,'qcat','qcat')]
    rho2 = results[(band,'ecat','qcat')]
    rho3 = results[(band,'wcat','wcat')]
    rho4 = results[(band,'qcat','wcat')]
    rho5 = results[(band,'ecat','wcat')]
    print('rho1.xip = ',rho1.xip)
    print('rho2.xip = ',rho2.xip)

    plt.clf()
    pretty_rho0(rho0.meanr, rho0.xip, np.sqrt(rho0.varxip))
    plt.savefig('rho0' + tag + '_' + band + '.pdf')

    plt.clf()
    pretty_rho1(rho1.meanr, rho1.xip, np.sqrt(rho1.varxip), 
                rho3.xip, np.sqrt(rho3.varxip),
                rho4.xip, np.sqrt(rho4.varxip))
    plt.savefig('rho1' + tag + '_' + band + '.pdf')

    plt.clf()
    pretty_rho2(rho2.meanr, rho2.xip, np.sqrt(rho2.varxip), 
                rho5.xip, np.sqrt(rho5.varxip))
    plt.savefig('rho2' + tag + '_' + band + '.pdf')
