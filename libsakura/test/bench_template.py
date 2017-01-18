#!/usr/bin/env python3
import timeit
from threading import Thread
from queue import Queue
import time

class MyBench:
	@staticmethod
	def createWS(id):
		return ()

	@staticmethod
	def destroyWS(id, ws):
		pass

	@staticmethod
	def createThreadLocal(thread_id):
		return (thread_id,)

	@staticmethod
	def destroyThreadLocal(tread_id, threadLocal):
		pass

	@staticmethod
	def run(tl, ws, id):
		time.sleep(0.01)

def worker(tl, context):
	while not context.stop:
		ws_id, ws = context.wss.get()
		if False:
			print("thr {0}: {1}".format(tl[0], ws_id))
		context.cls.run(tl, ws, ws_id)
		context.wss.put((ws_id, ws))
		context.count += 1

class Context:
	pass

def bench(cls, num_working_set, num_threads, numIters):
	assert num_working_set >= num_threads
	context = Context()
	context.cls = cls
	context.stop = False
	context.wss = Queue(num_working_set)
	for i in range(num_working_set):
		ws = cls.createWS(i)
		context.wss.put((i, ws))

	threads = []
	thread_locals = []
	context.count = 0
	for i in range(num_threads):
		tl = cls.createThreadLocal(i)
		thread_locals.append(tl)
		thr = Thread(target=worker, args=(tl, context))
		#thr.daemon = True
		threads.append(thr)

	start = time.monotonic()
	for thr in threads:
		thr.start()

	while context.count < numIters:
		time.sleep(1)
		now = time.monotonic()
		print("{0} workingSets/sec".format(context.count / (now - start)))

	context.stop = True

	for thr in threads:
		thr.join()

	for i, tl in enumerate(thread_locals):
		cls.destroyThreadLocal(i, tl)
	while not context.wss.empty():
		i, ws = context.wss.get()
		cls.destroyWS(i, ws)

if __name__ == "__main__":
	bench(MyBench, 4, 4, 2000)
