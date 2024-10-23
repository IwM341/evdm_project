from ._cpp_lib.pyevdm import *


def get_property(m_dict,address):
	if(len(address) == 0):
		return m_dict
	elif(len(address) == 1):
		return m_dict[address[0]]
	else:
		return get_property(m_dict[address[0]],address[1:])
def set_property(m_dict,address,value):
	if(address == 0):
		raise ValueError("in set_property expect len(address) > 0")
	elif(len(address) == 1):
		m_dict[address[0]] = value
	else:
		return set_property(m_dict[address[0]],address[1:],value)

evdm_restore_grid = {
	'evdm.Distrib':Distrib,
	'evdm.Capture':Capture,
	'evdm.Matrix':Matrix,
	'evdm.AnnPre':AnnPre}

evdm_restore= {
	'evdm.Distrib':scatter_event_info}

def _evdm_type_(m_object):
	return type(m_object) in [Distrib,Capture,Matrix,AnnPre,scatter_event_info]


def deep_move_impl(m_object,m_grids,m_bodies):
	if(type(m_object) == dict):
		if('type' in m_object):
			if(m_object['type'] in evdm_restore):
				return evdm_restore[m_object['type']](m_object)
			elif(m_object['type'] in evdm_restore_grid):
				return lambda mgrd: evdm_restore_grid[m_object['type']](mgrd,m_object)
			elif(m_object['type'] == '_evdm.BodyRef'):
				return m_bodies[m_object['index']]
			elif(m_object['type'] == '_evdm.GridRef'):
				return m_grids[m_object['index']][0]
		return {key : deep_move_impl(val,m_grids,m_bodies) for (key,val) in m_object.items()}
	if(type(m_object) == tuple or type(m_object) == list):
		return type(m_object)([deep_move_impl(val,m_grids,m_bodies) for val in m_object])
	else:
		return m_object

def deep_move(m_dict):
	m_bodies = []
	m_grids = []
	for bd in m_dict['bodies']:
		m_bodies.append( Body(bd) )
	for gr in m_dict['grids']:
		m_grids.append((GridEL(m_bodies[gr['bodyref']],gr['value']), gr['object_paths'] ))
	m_copy = deep_move_impl(m_dict,m_grids,m_bodies)
	for grd in m_grids:
		for m_path in grd[1]:
			m_constructor = get_property(m_copy,m_path)
			set_property(m_copy,m_path,m_constructor(grd[0]))
	return m_copy

def _add_body(_m_bodies,_body):
	_m_bodies.append(_body)
	return len(_m_bodies-1)
def _process_body(_m_bodies,_body):
	try:
		m_index = _m_bodies.inedex(_body)
		return m_index
	except:
		_m_bodies.append(_body)
		return len(_m_bodies)-1
def _process_grid(_m_grids,_m_bodies,_grid):
		index = 0
		for _m_g in _m_grids:
			if(_m_g['value'] == _grid):
				return index
			index += 1
		body_ref = _process_body(_m_bodies,_grid.body)
		_m_grids.append({'value':_grid,'bodyref':body_ref,'object_paths' : [] })
		return index


def _contain_evdm_(value):
	try:
		value._evdm_type_()
		return True
	except:
		pass
	
	try:
		m_props = value.__dict__
	except:
		raise TypeError("expect only classes in _contain_evdm_: error in serialize_Object")
	sum = False
	

def serialize_Object(m_object,m_grids,m_bodies,m_path = []):
	'''
		bodies = [body0,body1,...]
		grids = [{'value':grid0,'bodyref':0,object_paths = ["."] }]
	'''

	if(type(m_object) != dict):
		if(type(m_object) == Body):
			return {'type': '_evdm.BodyRef','index':_process_body(m_bodies,m_object)}
		elif(type(m_object) == GridEL):
			return {'type': '_evdm.GridRef','index':_process_grid(m_grids,m_bodies,m_object)}
		else:
			if(_evdm_type_(m_object)):
				if(hasattr(m_object,'grid')):
					v_grid = m_object.grid
					m_index = _process_grid(m_grids,m_bodies,v_grid)
					m_grids[m_index]['object_paths'].append(m_path)
					return m_object.to_object()
				else:
					return m_object.to_object()
			else:
				if(hasattr(m_object,'__dict__')):
					return serialize_Object(m_object.__dict__,m_grids,m_bodies,m_path)
				else:
					if(type(m_object) == tuple):
						return tuple([serialize_Object(m_object[i],m_grids,m_bodies,m_path + [i]) for i in range(len(m_object))])
					if(type(m_object) == list):
						return [serialize_Object(m_object[i],m_grids,m_bodies,m_path + [i]) for i in range(len(m_object))]
					if(type(m_object) == set):
						m_object_list = list(m_object)
						return [serialize_Object(m_object_list[i],m_grids,m_bodies,m_path + [i]) for i in range(len(m_object_list))]
					return m_object
	else:
		ret_dict = {}
		for (key,value) in m_object.items():
			n_path = m_path.copy()
			n_path.append(key)
			ret_dict[key] = serialize_Object(value,m_grids,m_bodies,n_path)
		return ret_dict