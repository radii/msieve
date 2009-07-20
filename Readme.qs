MSIEVE: A Library for Factoring Large Integers
Quadratic Sieve Module Documentation
Jason Papadopoulos


Introduction
------------

This file describes the implementation of the quadratic sieve that
Msieve uses. Whereas previous versions of Msieve could only use the
quadratic sieve, versions starting with v1.04 have the number field
sieve available, and that implementation is described in Readme.nfs

In factoring jargon, Msieve is an implementation of the self-initializing
multiple polynomial quadratic sieve (MPQS) with double large primes.
Unlike the NFS implementation, which is quite crude at the moment,
the QS implementation is very stable. To my knowledge it's the fastest
QS code available, sometimes by a huge margin. That's not to say that
it can't be better than it is now, but I've essentially run out of 
ideas to improve it. 

Since the QS implementation is quite good and the NFS implementation
is not, Msieve currently will always use the QS implementation by default.
It's early enough in the NFS development cycle that I have no precise
idea of the crossover point between QS and NFS in the specific case of 
the Msieve library; experience with other NFS implementations suggests
that it should eventually be slightly under the 100-digit problem size.
Nevertheless, the QS side is still useful at all problem sizes.

Since the library works essentially automatically once given a number
to factor, the only caveats on using the QS implementation involve 
cases where manual labor is required. To whit:


Distributed Computing
---------------------

By default, the demo application will keep going until all integers given to
it have been factored. If you happen to have several computers available and
need to factor a really big integer, it is possible to get all of the 
computers to work on factorizing the same number. This process is not automated
and is a little labor-intensive, so don't be surprised if it turns out to be
somewhat inconvienient. That's the price you pay for 1) keeping the code and
the user interface simple but 2) finishing your factorization much faster.

Before giving you recipes for 'going distributed', you have to know in general
terms how the quadratic sieve works. This is not the place to explain that
in detail (the references do it much better than I can) but you'll be a lot
less confused if you grasp the following facts:

	1. The quadratic sieve is a fast way of finding 'relations'.
	   A relation is a statement about the number you want to factor.

	2. A single relation is theoretically capable of factoring your
	   number. In practice that's impossible; the best you can do
	   is find a lot of different relations that fail. For a really 
	   big number, you'll need millions of these relations. The process
	   of finding relations is called 'sieving'.

	3. Once you have 'enough' of these relations, you can take the
	   millions of relations and combine them together into
	   a few dozen super-relations that *are* capable of finishing the
	   factorization with high probability.
	
	4. You don't know beforehand how many relations will be 'enough'. 
	   If you already have a collection of relations, you can use some 
	   math to determine whether there are 'enough' relations in that 
	   collection.

	5. The more sieving you do, the more relations you'll find. A block
	   of sieving work will produce an approximately constant number of
	   relations. So every computer you use for sieving will accumulate
	   relations at a constant rate. Obviously, the faster the computer
	   the higher the constant rate.

So (at least from low earth orbit) the quadratic sieve is divided into
two phases, the 'sieving' phase and the 'combining' phase. For msieve, the
combining phase takes a few minutes at most, but the sieving phase can take 
days or weeks. It's the sieving phase we want to do in parallel. We have to
keep sieving until all of the computers involved have produced a total of
'enough' unique relations.

Now that the preliminaries are out of the way, here are the recipes for
going distributed:

RECIPE 1:
   1. Start msieve on each sieving machine
   2. Every so often (say, once a day):
        - stop msieve on each sieving machine
	- combine the savefiles from all the sieving machines into
	  a single 'super-savefile'. Leave the savefiles from the
	  sieving machines alone otherwise
	- start one copy of msieve from the super-savefile. If it
	  finds there are 'enough' relations, then it will automatically
	  start the combining phase. 
	- If there are not 'enough' relations, the one copy of msieve
	  will start sieving again. Stop it, delete the super-savefile,
	  and repeat step 1.

RECIPE 2 (less data transfer):
   1. Start msieve on each sieving machine
   2. Every so often (say, once a day):
        - stop msieve on each sieving machine
	- APPEND the savefiles from all the sieving machines into
	  the existing 'super-savefile', then delete the savefiles 
	  from the sieving machines
	- start one copy of msieve from the super-savefile. If it
	  finds there are 'enough' relations, then it will automatically
	  start the combining phase. 
	- If there are not 'enough' relations, the one copy of msieve
	  will start sieving again. Stop it and repeat step 1. Do not
	  delete the super-savefile

Some notes on these recipes:

1. The copies of msieve that are sieving do not need to know about each other.
   Each one uses random numbers to begin the sieving process, and it's
   essentially impossible that two sieving machines will produce the same
   output by accident. Even if they did, the combining phase will notice
   and will remove the duplicates.

2. Each copy of msieve will give a running progress report of how much of the
   factorization it has completed. These reports assume there is only one
   sieving machine, so if you're using more than one sieving machine you 
   should not believe this output. When starting from the super-savefile,
   there really *is* one machine running so the output is accurate. If starting
   from the super-savefile shows that you only have a few relations left
   to go, you don't have to restart all the sieving machines; just let the one
   copy of Msieve finish the job.

3. At no point should the super-savefile contain any duplicate relations.
   Obviously, finding one relation and copying it a million times will not
   get you any closer to finishing the factorization. If you do this, you'll 
   fool msieve into thinking it has 'enough' relations, the combining phase
   will start, all the duplicates will get deleted and the combining phase
   will fail. It's just wasted effort. Likewise, sieving machines must all
   generate their own relations; DO NOT give relations back to sieving machines
   or share relations found across sieving machines. All you'll do is 
   generate a huge number of duplicates that will just get removed anyway.

4. Msieve uses the same name by default for its savefile. If your sieving
   machines write to a network directory, make sure they do not all write to
   the same file. There is no synchronization that keeps sieving machines
   from stomping on each other's output, and the code needed to lock the
   savefile is too machine-specific so I didn't add it. Msieve has a *lot*
   of error checking in the combining phase, so it's possible that the 
   factorization will succeed even if everybody did write to the same savefile. 
   However, every corrupted relation that must be skipped in the combining
   phase represents good work that you've wasted.

5. You will find for a really big factorization that relations start to
   accumulate grindingly slowly. This is normal; as more work gets done, the
   number of relations shoots up very quickly. By the time you're almost
   finished, relations will be collecting at incredible speed.

With the recipes in mind and a little patience, distributed sieving can be
a powerful tool for finishing your factorizations faster.
