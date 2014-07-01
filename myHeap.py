class MaxHeap:
	'''The goal of this was to be able to optimally remove zeros from a symmetric matrix while preserving 
	all of its intrinsic properties and maintaining symmetry.A list based binary max heap priority queue is used to 
	return the row/col to be removed. Its priority is determined based on the number of zeros it contains.
	 The heap is always left complete. The left child of a node
	at position p is at position 2p and the position of the right child is at 2p+1. I made my list structure 
	artificially 1-indexed to simplify integer division. A list-based hash table is used to keep track of where 
	a row/col is represented in the tree for constant time removal.'''
#------------------------------------------------------------------------------------#
	def __init__(self,matrix):
		self._dim=len(matrix)
		self._heapList=[float("inf")]
		self._hash = [None]*self._dim
		self._size=0
#------------------------------------------------------------------------------------#
	def fill(self,matrix):
		for i,row in enumerate(matrix):
			zeros=row.count(0)
			self.insert(zeros,i)

		assert None not in self._hash, "hashtable not filled properly"
		assert self._size==self._dim
#------------------------------------------------------------------------------------#
	def insert(self,zeros,index):
		'''insert something into the heap'''

		assert zeros>=0, "Please use non-zero input"
		assert index>=0, "Please use non-zero input"

		self._heapList.append((zeros,index))
		self._size += 1
		self._hash[index] = self._size
		self.upHeap(self._size)
#------------------------------------------------------------------------------------#
	def downHeap(self,position):
		'''restores heap order'''
		while position*2 <= self._size:
			maximum=self.maximum(position)
			if self._heapList[position]<self._heapList[maximum]:
				self.swap(position,maximum)
			position=maximum
#------------------------------------------------------------------------------------#
	def upHeap(self,position):
		'''restores heap order'''
		while position//2>0:
			if self._heapList[position][0]>self._heapList[position//2][0]:
				self.swap(position, position//2)
			position=position//2
#------------------------------------------------------------------------------------#
	def swap(self,index_a,index_b):
		'''swaps elements in the list'''

		hash_a=self._heapList[index_a][1]
		hash_b=self._heapList[index_b][1]

		old = self._heapList[index_a]
		self._heapList[index_a]=self._heapList[index_b]
		self._heapList[index_b]=old

		self._hash[hash_a]=index_b
		self._hash[hash_b]=index_a
#------------------------------------------------------------------------------------#
	def maximum(self,position):
		'''gets the min of a position's children'''
		if (position*2+1)>self._size or self._heapList[position*2]>self._heapList[position*2+1]:
				return position*2
		return (position*2+1)
#------------------------------------------------------------------------------------#
	def remove_max(self):
		'''removes and returns the maximum of the heap. I'm being lazy and using assertions instead
		of throwing my own errors but it will throw an assertion error if the heap is empty.'''
		assert self._size>1, "heap is empty"
		heap_max = self._heapList[1]
		index=heap_max[1]
		self.swap(1,self_size)
		self._size-=1

		for heap_pos in self._hash[index+1:]:
			self._heapList[heap_pos]=(self._heapList[heap_pos][0],self._heapList[heap_pos][1]-1)
			assert self._heapList[heap_pos][1]>0

		
		self._hash.pop(index)
		self._heapList.pop()
		self.downHeap(1)
		return heap_max

#------------------------------------------------------------------------------------#

	def remove_zero(self,index):
		heap_position = self._hash[index]
		num_zeros, i = self._heapList[heap_position]
		self._heapList[heap_position]=(num_zeros-1,i)
		self.downHeap(heap_position)

