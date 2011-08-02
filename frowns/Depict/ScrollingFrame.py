from Tkinter import *
from Tkinter import _cnfmerge

# This object acts very much like a Frame, but if it gets too big to display on
# the screen, it creates scrollbars that let the user scroll around it's viewing
# area.
class ScrollingFrame(Frame):
    def __init__(self, parent=None, **kw):
        # Create a canvas for this frame to appear in.  The canvas is contained in
        # the ScrollingFrame object, but the TK frame will appear in the canvas.
        self.container = Frame(parent)
        self.canvas = Canvas(self.container)
        self.scrollX = None
        self.scrollY = None

        # Now initialize the frame with the container as it's parent (sneaky...)
        apply(Frame.__init__, (self, self.container), kw)

        # Allow frame resizes to re-init scrolling frame
        self.container.bind("<Configure>", self.ResizeWindowFunc, "+")

    def ResizeWindowFunc(self, event):
        try:
            self.update()
        except:
            pass

    # Redefine update() to create scrollbars if frame size exceeds window size
    def update(self):
        # Call parent class update function for ScrollingFrame
        Frame.update(self)
        self.container.update()

        frame_width = self.winfo_width()
        frame_height = self.winfo_height()

        # After updating frame (and storing dimension values), configure canvas to
        # frame size (This sets up proper default dimensions, to fit frame snugly)
        # Only do this if size doesn't already agree
        c_width = self.canvas.winfo_width()
        c_height = self.canvas.winfo_height()
        if c_width != frame_width:
            self.canvas.configure(width = frame_width)
        if c_height != frame_height:
            self.canvas.configure(height = frame_height)

        # Then pack and update canvas. Container frame will be clipped to screen
        # dimensions, but the Canvas and ScrollingFrame will not.
        # This allow us to test if scrollbars should be displayed, by seeing if
        # ScrollingFrame is bigger than container frame
        self.canvas.grid(column=1, row=1, sticky=N+S+E+W)
        self.canvas.update()

        container_width = self.container.winfo_width()
        container_height = self.container.winfo_height()

        # Use the narrower of the two, frame or container dimension, as the new
        # container dimension (This is because scrollbars are in the container, and
        # add to it's size; which looks ugly when canvas is constantly resized between
        # alternating values).  Failure to do this step can cause some wierd
        # oscillating as canvas and container fight to resize themselves correctly
        container_width = min(container_width, frame_width)
        container_height = min(container_height, frame_height)

        # Now test to see if scrollbars need to be created.  If either has happened,
        # container dimensions have changed and must be reset
        if frame_width > container_width and not self.scrollX:
            # Create horizontal scrollbar and tie it to canvas
            self.scrollX = Scrollbar(self.container, orient=HORIZONTAL)
            self.canvas['xscrollcommand'] = self.scrollX.set
            self.scrollX['command'] = self.canvas.xview
            self.scrollX.grid(column=1, row=0, sticky=E+W)
        # Otherwise, test to see if scrollbars need to be destroyed
        elif frame_width <= container_width and self.scrollX:
            self.scrollX.grid_forget()
            self.scrollX = None

        if frame_height > container_height and not self.scrollY:
            # Create vertical scrollbar and tie it to canvas
            self.scrollY = Scrollbar(self.container, orient=VERTICAL)
            self.canvas['yscrollcommand'] = self.scrollY.set
            self.scrollY['command'] = self.canvas.yview
            self.scrollY.grid(column=0, row=1, sticky=N+S)
        elif frame_height <= container_height and self.scrollY:
            self.scrollY.grid_forget()
            self.scrollY = None

        # If we are using scrollbars, set canvas size to window dimensions, and
        # scrollregion to ScrollingFrame size.  Now scrollbars should work fine
        if self.scrollX or self.scrollY:
            self.canvas.configure(width = container_width, height = container_height)
            self.canvas.configure(scrollregion = (0, 0, frame_width, frame_height))

    # Force pack to pack both canvas and frame
    def pack(self, cnf = {}, **kw):
        if kw:
            cnf = _cnfmerge((cnf, kw))

        # Add frame to the canvas AFTER packing container frame
        # (to set width and height properly)
        self.container.pack(cnf)
        try:
            self.canvas.delete(self.canvas_frame_id)
        except AttributeError:
            pass
        self.canvas_frame_id = self.canvas.create_window(0, 0, window=self, anchor=NW)

    # Force bind to bind to container frame, not internal frame (which would probably
    # cause a crash)
    def bind(self, cnf = {}, **kw):
        if kw:
            cnf = _cnfmerge((cnf, kw))
        self.container.bind(cnf)

    # Force bind to bind to container frame, not internal frame (which would probably
    # cause a crash)
    def unbind(self, cnf = {}, **kw):
        if kw:
            cnf = _cnfmerge((cnf, kw))
        self.container.unbind(cnf)

