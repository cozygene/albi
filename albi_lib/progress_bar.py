# Progress bar

import sys, os

SHELL_WIDTH = 70

def ioctl_GWINSZ(fd):
    try:
        import fcntl, termios, struct
        _, width = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,'1234'))
    except:
        return
    return width

def get_shell_width():
    width = ioctl_GWINSZ(0)

    if not width:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            width = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass

    if not width:
        width = os.environ.get('COLUMNS')

    return width if width else SHELL_WIDTH


class ProgressBarIter(object):
    def __init__(self, length, stdout = sys.stdout, fill = '#', width = get_shell_width()):
        self.length = float(length)
        if stdout.isatty():
          self.stdout = stdout
        else:
          self.stdout = sys.stdout
        self.current = 0
        self.fill = fill
        self.width = width
        self.prefix = '| '
        self.suffix = ' |'
        self.percentage = " {percentage:.2f}%"

    def __iter__(self):
        return self

    def __next__(self): # Python 3 uses __next__ for iterator
        return self.next()

    def next(self): # Python 2 uses next for iterator
        percentage = (self.current / self.length)
        percentage_str = self.percentage.format(percentage = percentage * 100)
        width_left = self.width - len(self.prefix) - len(self.suffix) - len(percentage_str)
        process_bar_fill_str = ('#' * int(percentage * width_left)).ljust(width_left)
        self._write(self.prefix + process_bar_fill_str + percentage_str + self.suffix)

        if percentage == 1:
            raise StopIteration

        self.current += 1
        return percentage

    def _write(self, line) :
        self.stdout.write('\r') # cleans up the stdout
        self.stdout.write(line)
        self.stdout.flush()

    def __len__(self):
      return int(self.length)
