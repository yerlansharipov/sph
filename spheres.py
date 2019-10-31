#!/anaconda3/bin/python

# Copyright 2019 yerlan sharipov@bu.edu
import numpy as np
import math
import sys
import time
import multiprocessing
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Spheres():

	def __init__ (self, name='newbie', mass=0, rad=0, pos=(0,0,0), vel=(0,0,0), col_times = -1, cont=False, shock=False):
		self.name = name
		self.rad = rad
		self.mass = mass
		self.pos = pos
		self.vel = vel
		self.cont = cont
		self.col_times = col_times
		self.shock = shock

	def collision_detection(self, other):

		col_time = -1
		col_pos = np.dot(np.subtract(self.pos, other.pos),np.subtract(self.vel, other.vel))

		if col_pos < 0 or other.cont == True or self.cont == True:
			vi_vj = np.subtract(self.vel, other.vel)
			pi_pj = np.subtract(self.pos, other.pos)
			vi_vj_sq = np.dot(vi_vj, vi_vj)
			pi_pj_sq = np.dot(pi_pj, pi_pj)
			a = 0
			b = 0
			c = 0

			a = vi_vj_sq
			b = 2 * np.dot(vi_vj, pi_pj)
			if other.cont == True or self.cont == True:
				c = pi_pj_sq - (self.rad - other.rad)**2
			else:
				c = pi_pj_sq - (self.rad + other.rad)**2

			discr = b**2 - 4*a*c
			if discr > 0:
				# col_time = min(abs((-b - math.sqrt(discr)) / (2*a)), abs((-b + math.sqrt(discr)) / (2*a)))
				ct_1 = ((-b - math.sqrt(discr)) / (2*a))
				ct_2 = ((-b + math.sqrt(discr)) / (2*a))
				if ct_1 > 0 and ct_2 > 0:
					col_time = min (ct_1, ct_2)
				elif ct_1 > 0:
					col_time = ct_1
				elif ct_2 > 0 :
					col_time = ct_2

				if round(col_time,2) == 0:
					col_time = max(abs((-b - math.sqrt(discr)) / (2*a)), abs((-b + math.sqrt(discr)) / (2*a)))
			elif discr < 0:
				col_time = -1
			elif discr == 0:
				if a != 0:
					col_time = (-b) / (2*a)
				else:
					col_time = -1
		else:
			discr = -1
			col_time = -1
			a= 0
			b=0
			c=0
		print("[DEBUG] NAMES", self.name, other.name, "|| COL TIME: ", col_time, "ABC:", a,b,c, " || discr: ", discr)
		return col_time

	def ellastic_collision(self, other, time_step):
		
		# updates velocities after collision

		if self.cont == False and other.cont == False:
			r1_r2 = np.subtract(self.pos, other.pos)
			v1_v2 = np.subtract(self.vel, other.vel)
			r1_r2_sq = np.dot(r1_r2,r1_r2)
			mms = self.mass + other.mass
			denom = np.dot(mms, r1_r2_sq)
			nom = np.dot(v1_v2, r1_r2)
			nom = np.dot(nom, r1_r2)
			nom = np.dot(nom, 2)
			self.vel = (np.subtract(self.vel, other.mass*nom/denom))
			other.vel = (np.add(other.vel, self.mass*nom/denom))
		else:

			# delta_pos = sqrt( other.rad**2 / (np.dot(other.pos[0], other.pos[0])  + np.dot(other.pos[1], other.pos[1]) + np.dot(other.pos[0], other.pos[0])))

			# x = np.dot(delta_pos,other.pos[0])
			# y = np.dot(delta_pos,other.pos[1])
			# z = np.dot(delta_pos,other.pos[2])

			hit_pos = np.dot(other.pos, self.rad/(self.rad-other.rad))

			print("POSITON: ", other.pos, "DOT POSITION: ", hit_pos)
			r1_r2 = -np.subtract(tuple(hit_pos), other.pos)
			v1_v2 = -np.subtract(self.vel, other.vel)
			r1_r2_sq = np.dot(r1_r2,r1_r2)
			mms = 1
			denom = np.dot(mms, r1_r2_sq)
			nom = np.dot(v1_v2, r1_r2)
			nom = np.dot(nom, r1_r2)
			nom = np.dot(nom, 2)
			#self.vel = (np.subtract(self.vel, other.mass*nom/denom))
			other.vel = (np.subtract(other.vel, nom/denom))


			# other.vel = np.dot(-1,other.vel)


		# set shock flag after collision
		if self.shock == False and self.name != 'The Universe':
			self.shock = True
		elif other.shock == False and other.name != 'The Universe':
			other.shock = True
			
		#decrement col_times
		self.col_times -= 1
		other.col_times -= 1


################################################################################################################
################################################################################################################


def process_input():

	print("Please enter the mass, radius, x/y/z position, x/y/z velocity and name of each sphere. When complete, use EOF / Ctrl-D to stop entering")

	big_rad = int(sys.argv[1])
	col_times = int(sys.argv[2])

	spheres_parameters = [x.strip() for x in sys.stdin.readlines()]
	
	spheres = dict()
	energy = 0
	momentum = 0

	# initialize The Universe and spheres
	print("Here are the initial spheres.")
	spheres['The Universe'] = Spheres(name='The Universe', rad = big_rad, cont=True)
	for b in spheres_parameters:
		b = b.split()
		spheres[b[8]] = Spheres(name=b[8], rad = int(b[1]), mass=int(b[0]), pos=(int(b[2]), int(b[3]), int(b[4])), vel=(int(b[5]), int(b[6]), int(b[7])), col_times = col_times)

		energy += 0.5 * spheres[b[8]].mass * np.dot(spheres[b[8]].vel, spheres[b[8]].vel)
		momentum += np.dot(spheres[b[8]].mass, spheres[b[8]].vel)

		print("%-10s m=%-4s R=%-4s p=%-20s v=%-10s" % (spheres[b[8]].name, spheres[b[8]].mass, spheres[b[8]].rad, spheres[b[8]].pos, spheres[b[8]].vel))

	print("energy:   %-5s" % (energy))
	print("momentum: %-10s" % (str(tuple(momentum))))

	return spheres

def next_collision(spheres):

	# calculate collision time for each sphere
	collisions = {}
	spheres_list = list(spheres)

	for key in list(spheres):
		spheres_list.pop(0)
		for key2 in spheres_list:
			col_time = spheres[key].collision_detection(spheres[key2])
			# print("shock: ", spheres[key].name, spheres[key].shock, spheres[key2].name, spheres[key2].shock)
			if col_time > 0:
				collisions[(key , key2)] = col_time
	# print("[INFO] COLLISIONS: ", collisions)

	# check the closest collisions
	closest_col_time = min(collisions.values())
	colliding_objects = []
	colliding_objects = [ key for key, val in collisions.items() if val == closest_col_time ]

	return colliding_objects, closest_col_time

def plots(spheres):

	#plot 2-dimensional
	fig, ax = plt.subplots()

	for sphere in spheres.items():
		if sphere[1].name == "The Universe":
			ax.add_patch(plt.Circle((sphere[1].pos[0], sphere[1].pos[1]), sphere[1].rad, color='r', alpha=0.1))
			ax.set_aspect('equal', adjustable='datalim')
		else:	
			ax.add_patch(plt.Circle((sphere[1].pos[0], sphere[1].pos[1]), sphere[1].rad, color='r', alpha=0.5))
			ax.set_aspect('equal', adjustable='datalim')

	ax.plot()
	plt.show()

def plot3d(spheres):
	fig = plt.figure()
	ax = plt.axes(projection="3d")
	for sphere in spheres.items():
		if sphere[1].name != "The Universe":
			ax.scatter3D(sphere[1].pos[0], sphere[1].pos[1])
		

	plt.show()

def main():

	# process input from terminal and return dictionary of spheres {name: object}
	spheres = process_input()
	print("=====================================================================================================")
	print("\n")	


	# plot3d(spheres)
	plots(spheres)
	# create The Universe and maintain it untill objects dissapear
	while len(spheres) > 1:

		# check for the next collisions and return elements' pairs which are going to collide and time
		colliding_objects, closest_col_time = next_collision(spheres)
		
		# update position of all objects 
		print("[INFO]",colliding_objects, " collided after",round(closest_col_time,1), "s: recalculating system's parameters")
		for name, sphere in spheres.items():
			sphere.pos = tuple(np.add(sphere.pos, np.dot(closest_col_time, sphere.vel)))
		
		# perform elastic collision of objects
		for objects in colliding_objects:
			obj_1 = spheres[objects[0]]
			obj_2 = spheres[objects[1]]
			obj_1.ellastic_collision(obj_2, closest_col_time)

		# show debug info
		for name, sphere in spheres.items():
			if name != "The Universe":
				print("%-15s m=%-4s R=%-4s p=%-20s v=%-10s" % (sphere.name, sphere.mass, sphere.rad, (round(sphere.pos[0],1), round(sphere.pos[1],1), round(sphere.pos[2],1)), (round(sphere.vel[0],1), round(sphere.vel[1],1), round(sphere.vel[2],1))))

		# call plot method
		plots(spheres)

		# delete objects with expired number of collisions
		for name in list(spheres):
			if spheres[name].col_times == 0:
				del spheres[name]
				print("[INFO]", name, "is gonna be deleted")

		print()
		print("=====================================")
		print()
		
if __name__ == '__main__':
	main()