
import multiprocessing
import segfault
import signal
import time

if __name__ == '__main__':

    def this_dies():
        def sig_handler(signum, frame):
            print('segfault')

        signal.signal(signal.SIGSEGV, sig_handler)

        print('Hello World')
        time.sleep(1)
        segfault.segfault()

    p = multiprocessing.Process(target=this_dies())
    p.start()
    try:
        p.join()
    except Exception as e:
        pass

    print('I am joined')
