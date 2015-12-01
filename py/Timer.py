__author__ = 'Jwely'

from datetime import datetime



class Timer:

    def __init__(self):
        self.start_time = datetime.now()
        self.laps = []
        self.finish_time = None

    def lap(self):
        lap = datetime.now()
        self.laps.append(lap)
        delta = lap - self.start_time
        return delta.seconds + delta.microseconds / 1.0e6

    def finish(self):
        self.finish_time = datetime.now()
        delta = self.finish_time - self.start_time
        return delta.seconds + delta.microseconds / 1.0e6


if __name__ == "__main__":
    import time

    t = Timer()
    time.sleep(1.128)
    print t.finish()
