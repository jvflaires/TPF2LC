"""
Code to create a mask and plot the corresponding LC from a 
lightkurve Target Pixel File (TPF) object
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from matplotlib.gridspec import GridSpec
from astropy.time import Time
import lightkurve as lk



# Main Class
class TPF2LC:
    """
    Class to interactively select the star pixels 
    in a Lightkurve Target Pixel File (TPF).
    """

    # Class initialization
    def __init__(self, tpf):
        """
        Initialize the Class object with a Target Pixel File (TPF).

        Parameters:
        -----------
        tpf : lightkurve.TargetPixelFile
            Object's Target Pixel File (TPF) obtained from lightkurve package.
        """
        self.tpf = tpf  # Target Pixel File          
        self.data = np.nansum(self.tpf.flux, axis=0) / np.median(np.nansum(self.tpf.flux, axis=0)) # 2D flux normalized data
        self.current_x, self.current_y = 0, 0  # Pixel coordinates
        self.marked_pixels = []  # List of marked pixels
        self.hashed_pixels = []  # list of hashed squares
        self.selected_pixels = np.zeros(self.data.shape, dtype='bool')  # Matrix of selected pixels
        self.fig, self.ax1, self.ax2 = None, None, None  # Matplotlib figure and axis
        self.rect = None  # Hashed square
        self.text = None  # Show pixel values


    # Function that display the build-in functions of the Class
    def help(self):
        """
        Description of the build-in functions of the Class
        """
        print("Name: TPF2LC")
        print("Description: Interactively creates a Target Pixel File (TPF) mask and")
        print("             produces the corresponding light curve.                 ")
        print("Author: João Aires (UFRN)                                            ")
        print("                                                                     ")
        print("=====================================================================")
        print("Function                      Return                                 ")
        print("=====================================================================")
        print("interact()                    Display the interactive module to select")
        print("                              the mask.                              ")
        print("get_mask()                    Returns the mask created in interact().")
        print("get_lightcurve()              Returns the light curve assosciated with")
        print("                              the mask selection.                    ")
        print("get_flatten_lightcurve()      Returns the flattened light curve assos-")
        print("                              ciated with the mask selection.        ")
        print("=====================================================================")
        print("                                                                     ")
        print("IMPORTANT: You need to use %matplotlib QtAgg to interaction work.    ")

    
    # Function to update square position in pixels
    def update_rectangle(self):
        """
        Update the rectangle that runs through the image pixels.
        """
        self.rect.set_xy((self.current_x - 0.5, self.current_y - 0.5))                         
        self.text.set_text(f"Pixel ({self.current_x}, {self.current_y}): {self.data[self.current_y, self.current_x]:.2f}")
        self.text.set_position((self.current_x, self.current_y - 1))                           
        self.fig.canvas.draw()


    # Function to update hashed squares
    def update_hashed_pixels(self):
        """
        Update the hashed pixels selected.
        """
        for mark in self.hashed_pixels:
            mark.remove()
        self.hashed_pixels.clear()
        for (x, y) in self.marked_pixels:
            mark = plt.Rectangle((x-0.5, y-0.5), 1, 1, edgecolor='white', facecolor='None', hatch='//', alpha=1)
            self.ax1.add_patch(mark)
            self.hashed_pixels.append(mark)
        self.fig.canvas.draw()


    # Function to update the light curve
    def update_lc(self):
        """
        Updates the light curve during the pixel selection to the one produced
        by thoses pixels that were selected.
        """
        self.ax2.clear()

        # Check if the is any selected pixels to plot the corresponding Light curve
        if self.selected_pixels.any() == True:
            lc = self.tpf.to_lightcurve(aperture_mask=self.selected_pixels)
            time = Time(lc.time)
            time_JD = time.jd
            flux = lc.flux.value
            flux_err = lc.flux_err.value
            self.ax2.errorbar(time_JD-2450000, flux / np.nanmedian(flux), yerr=flux_err/np.nanmedian(flux),
                              marker='.', ls='', color='k', lw=.7, capsize=1)
            self.ax2.set_xlabel('Time (JD)')
            self.ax2.set_ylabel('Normalized flux')
            self.fig.canvas.draw()
        else:
            self.ax2.plot()
            self.ax2.set_xlabel('Time (JD)')
            self.ax2.set_ylabel('Normalized flux')
            self.fig.canvas.draw()


    # Function defining the events
    def on_key(self, event):
        """
        Reads the keyboard events and execute the corresponding actions
        """
        # Commands to run through the pixels
        if event.key == 'right':  # Move to the right
            self.current_x = min(self.current_x + 1, self.data.shape[1] - 1)
        elif event.key == 'left':  # Move to the left
            self.current_x = max(self.current_x - 1, 0)
        elif event.key == 'up':  # Move up
            self.current_y = min(self.current_y + 1, self.data.shape[0] - 1)
        elif event.key == 'down':  # Move down
            self.current_y = max(self.current_y - 1, 0)

        # Select the pixel
        elif event.key == '1':
            if (self.current_x, self.current_y) not in self.marked_pixels:
                self.marked_pixels.append((self.current_x, self.current_y))
                self.update_hashed_pixels()
                self.selected_pixels = np.zeros(self.data.shape, dtype='bool')
                for (x, y) in self.marked_pixels:
                    self.selected_pixels[y,x] = True
                self.update_lc()
                
        # Unselect the pixel
        elif event.key == '2':
            if (self.current_x, self.current_y) in self.marked_pixels:
                self.marked_pixels.remove((self.current_x, self.current_y))
                self.update_hashed_pixels()
                self.selected_pixels = np.zeros(self.data.shape, dtype='bool')
                for (x, y) in self.marked_pixels:
                    self.selected_pixels[y,x] = True
                self.update_lc()
            else:
                pass
                
        # End the pixel selection
        elif event.key == 'q':  # Sair
            plt.close()
            for (x, y) in self.marked_pixels:
                self.selected_pixels[y,x] = True
        
        self.update_rectangle()


    # Function that obtains the mask produced
    def get_mask(self):
        """
        Returns the final mask produced during selection.

        Returns:
        --------
        mask : numpy.ndarray
            Boolean mask where 'True' represents those pixels selected.
        """
        mask = np.zeros(self.data.shape, dtype='bool')
        for (x,y) in self.marked_pixels:
            mask[y,x] = True
        return mask


    # Function that returns the Light Curve obtained
    def get_lightcurve(self):
        """
        Returns the normalized light curve produced by applying the mask on the TPF.

        Returns:
        --------
        lc : lightkurve.LightCurve
            Light curve produced by the mask generated during selection.
        """
        mask = self.get_mask()
        lc = self.tpf.to_lightcurve(aperture_mask=mask).normalize()
        lc = lc['time', 'flux', 'flux_err']
        lc['time']     = lc['time'].jd
        lc['flux']     = lc['flux'].value
        lc['flux_err'] = lc['flux_err'].value
        return lc


    # Function to produce flatten light curve
    def get_flatten_lightcurve(self, window_length=101):
        """
        Returns the flatten and normalized light curve produced by applying the mask on the TPF.

        Returns:
        --------
        flat_lc : lightkurve.LightCurve
            Flatten light curve produced by the mask generated during selection.
        """
        mask = self.get_mask()
        flat_lc = self.tpf.to_lightcurve(aperture_mask=mask).normalize().flatten(window_length=window_length)
        flat_lc = flat_lc['time', 'flux', 'flux_err']
        flat_lc['time']     = flat_lc['time'].jd
        flat_lc['flux']     = flat_lc['flux'].value
        flat_lc['flux_err'] = flat_lc['flux_err'].value
        return flat_lc
        
    
    # Function to create a aperture mask
    def interact(self):
        """
        Starts the interaction with the TPF.
        """
        # Create figure object
        self.fig = plt.figure(figsize=(8,8))
        gs = GridSpec(3, 1, figure=self.fig)

        # Target Pixel File axis
        self.ax1 = self.fig.add_subplot(gs[0:2, 0])
        im1 = self.ax1.imshow(self.data, cmap='viridis', origin='lower')
        plt.colorbar(im1, ax=self.ax1, label=r"Flux ($\rm e^{-}$ / s)")
        self.ax1.set_title("'Arrows': navegate;     'q': close; \n '1': add mask,     '2': remove mask")  
        self.rect = plt.Rectangle((-0.5, -0.5), 1, 1, edgecolor='red', facecolor='none')
        self.text = self.ax1.text(0.5, -0.5, "", fontsize=12, color='red')
        self.ax1.add_patch(self.rect)

        # Light curve axis
        self.ax2 = self.fig.add_subplot(gs[2,0])
        self.ax2.set_xlabel('Time (JD)')
        self.ax2.set_ylabel('Normalized flux')
        
        # Connect to allow key event
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        
        # Show plot
        plt.show()