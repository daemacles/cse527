import cse527
import time
import sys
import heapq
import zmq
from multiprocessing import Process
import numpy as np
import base64
import util

import pyximport; pyximport.install()
from inner_loops import seqRC

accum_addr = "tcp://diglett.cs.washington.edu"
accum_port = '5598'
control_addr = "tcp://diglett.cs.washington.edu"
control_port = '5599'

def worker(num_bp):
    genome = cse527.Genome('sacCer3_yeast2011.2bit')
    context = zmq.Context()

    controller = context.socket(zmq.SUB)
    controller.connect(control_addr+':'+control_port)
    controller.setsockopt(zmq.SUBSCRIBE, "")
    
    accumulator = context.socket(zmq.PUSH)
    accumulator.connect(accum_addr+':'+accum_port)

    running = True
    while True:
        try:
            ctrl_message = controller.recv(zmq.NOBLOCK)
            if ctrl_message == 'STOP':
                context.destroy(linger=0)
                break
        except zmq.core.error.ZMQError, e:
            pass
        
        guess_hits, guess = cse527.runGuess(genome, num_bp)
        guess = base64.b64encode(guess.tostring())
        answer = {'guess_hits': guess_hits, 'guess':guess}
        accumulator.send_json(answer)
    return
        
def accumulator(iterations, num_guesses):
    best_guesses = [(0, 0, np.array([0], dtype=np.uint8))
                    for i in range(num_guesses)]
    heapq.heapify(best_guesses)
    insertGuess = heapq.heappushpop
    uid = 1       # unique id used to disambiguate identically scoring guesses

    context = zmq.Context()
    
    accumulator = context.socket(zmq.PULL)
    accumulator.bind("tcp://*:" + accum_port)

    controller = context.socket(zmq.PUB)
    controller.bind("tcp://*:" + control_port)

    step = iterations / 40.0
    next_level = step
    print 'Progress [ 0                                    ] %d'%iterations
    print '         ',
    for task_num in range(iterations):
        if task_num > next_level:
            sys.stdout.write('^')
            sys.stdout.flush()
            next_level += step
        ans_dict = accumulator.recv_json()
        guess = np.fromstring(base64.b64decode(ans_dict['guess']),
                              dtype=np.uint8)
        insertGuess(best_guesses, (ans_dict['guess_hits'], uid, guess))
        uid += 1
    print
        
    controller.send('STOP')
    accumulator.close()

    return best_guesses

def runServer(iterations, num_saved):
        start_time = time.time()
        best_guesses = accumulator(iterations, num_saved)
        time_taken = time.time() - start_time
        print 'That took %.2f seconds'%(time_taken)
        print '=========== Best guesses ============'
        print ' Hits           Gene+           Gene-'
        for hits, uid, guess in best_guesses:
            print '%5d %15s %15s'%(
                hits,
                util.seqToBaseString(guess),
                util.seqToBaseString(seqRC(guess)))
        return best_guesses

if __name__ == '__main__':
    def usage():
        print 'Usage:'
        print '  %s w <Number of worker processes> <num basepairs in guess>'%(sys.argv[0])
        print '-or-'
        print '  %s s <Number of guesses> <Number of guesses to save>'%(sys.argv[0])
        print 'At least one instance of a server ("s" option) must be launched.'
        return
    
    try:
        if sys.argv[1] == 'w':
            num_workers = int(sys.argv[2])
            num_bp = int(sys.argv[3])
            workers = [Process(target=worker, args=(num_bp,))
                       for work_num in range(num_workers)]
            for worker in workers: worker.start()

        elif sys.argv[1] == 's':
            iterations = int(sys.argv[2])
            num_saved = int(sys.argv[3])
            best_guesses = runServer(iterations, num_saved)
        else:
            usage()
            sys.exit(1)
    except IndexError:
        usage()
        sys.exit(1)
    except ValueError:
        usage()
        sys.exit(1)


