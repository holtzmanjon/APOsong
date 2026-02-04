from holtztools import plots 

def plot(wav,waves) :
    fig,ax=plots.multi(1,2,hspace=0.001)
    pix=wav.pix[wav.weights>0].astype(int)
    y=wav.y[wav.weights>0]
    d=waves[y,pix]-waves[y,pix+1]
    ax[0].scatter(wav.waves[wav.weights>0],wav.fwhm[wav.weights>0],c=wav.pix[wav.weights>0],s=5)
    im=ax[1].scatter(wav.waves[wav.weights>0],wav.waves[wav.weights>0]/(d*wav.fwhm[wav.weights>0]),c=wav.pix[wav.weights>0],s=5)
    ax[1].set_xlabel('Wavelength')
    ax[1].set_ylabel('Resolution')
    ax[0].set_ylabel('FWHM (2x2 pixels)')
    ax[0].set_ylim(0,10)
    ax[1].set_ylim(0,120000)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    cbar_ax.set_ylabel('X pixel')

