import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Flatten, Add,Lambda,LeakyReLU
from tensorflow.keras.layers import Dropout, Conv1D, MaxPooling1D,MaxPooling1D,GlobalMaxPooling1D,SpatialDropout1D
from tensorflow.keras.layers import Activation,concatenate,BatchNormalization
from tensorflow.keras.optimizers import SGD,Adam,schedules
from tensorflow.keras.utils import plot_model


from tensorflow.keras import regularizers
from tensorflow.keras import initializers
from tensorflow.keras import constraints
from tensorflow.keras import backend as K
#from tensorflow.keras.engine import Layer, InputSpec
from tensorflow.keras.layers import Layer, InputSpec
from tensorflow.keras.metrics import binary_accuracy
from tensorflow.keras.initializers import Ones, Zeros
#from tensorflow.keras.utils.generic_utils import get_custom_objects

# from tensorflow.keras.models import Sequential, Model
# from tensorflow.keras.layers.convolutional import MaxPooling1D, Convolution1D
# from tensorflow.keras.layers.core import Dense, Dropout, Activation, Flatten



class GroupNormalization(Layer):
	def __init__(self,groups=32,axis=-1,epsilon=1e-5,center=True,scale=True,beta_initializer='zeros',gamma_initializer='ones',beta_regularizer=None,
				 gamma_regularizer=None,beta_constraint=None,gamma_constraint=None,**kwargs):
		super(GroupNormalization, self).__init__(**kwargs)
		self.supports_masking = True
		self.groups = groups
		self.axis = axis
		self.epsilon = epsilon
		self.center = center
		self.scale = scale
		self.beta_initializer = initializers.get(beta_initializer)
		self.gamma_initializer = initializers.get(gamma_initializer)
		self.beta_regularizer = regularizers.get(beta_regularizer)
		self.gamma_regularizer = regularizers.get(gamma_regularizer)
		self.beta_constraint = constraints.get(beta_constraint)
		self.gamma_constraint = constraints.get(gamma_constraint)

	def build(self, input_shape):
		dim = input_shape[self.axis]

		if dim is None:
			raise ValueError('Axis '+str(self.axis)+' of input tensor should have a defined dimension but the layer received an input with shape '+str(input_shape)+'.')

		if dim < self.groups:
			raise ValueError('Number of groups ('+str(self.groups)+') cannot be more than the number of channels ('+str(dim)+').')

		if dim % self.groups != 0:
			raise ValueError('Number of groups ('+str(self.groups)+') must be a multiple of the number of channels ('+str(dim)+').')

		self.input_spec = InputSpec(ndim=len(input_shape),axes={self.axis: dim})
		shape = (dim,)

		if self.scale:
			self.gamma = self.add_weight(shape=shape,name='gamma',initializer=self.gamma_initializer,regularizer=self.gamma_regularizer,constraint=self.gamma_constraint)
		else:
			self.gamma = None

		if self.center:
			self.beta = self.add_weight(shape=shape,name='beta',initializer=self.beta_initializer,regularizer=self.beta_regularizer,constraint=self.beta_constraint)
		else:
			self.beta = None

		self.built = True

	def call(self, inputs, **kwargs):
		input_shape = K.int_shape(inputs)
		tensor_input_shape = K.shape(inputs)

		# Prepare broadcasting shape.
		reduction_axes = list(range(len(input_shape)))
		del reduction_axes[self.axis]
		broadcast_shape = [1] * len(input_shape)
		broadcast_shape[self.axis] = input_shape[self.axis] // self.groups
		broadcast_shape.insert(1, self.groups)

		reshape_group_shape = K.shape(inputs)
		group_axes = [reshape_group_shape[i] for i in range(len(input_shape))]
		group_axes[self.axis] = input_shape[self.axis] // self.groups
		group_axes.insert(1, self.groups)

		# reshape inputs to new group shape
		group_shape = [group_axes[0], self.groups] + group_axes[2:]
		group_shape = K.stack(group_shape)
		inputs = K.reshape(inputs, group_shape)

		group_reduction_axes = list(range(len(group_axes)))
		group_reduction_axes = group_reduction_axes[2:]

		mean = K.mean(inputs, axis=group_reduction_axes, keepdims=True)
		variance = K.var(inputs, axis=group_reduction_axes, keepdims=True)
		inputs = (inputs - mean) / (K.sqrt(variance + self.epsilon))

		# prepare broadcast shape
		inputs = K.reshape(inputs, group_shape)
		outputs = inputs

		# In this case we must explicitly broadcast all parameters.
		if self.scale:
			broadcast_gamma = K.reshape(self.gamma, broadcast_shape)
			outputs = outputs * broadcast_gamma

		if self.center:
			broadcast_beta = K.reshape(self.beta, broadcast_shape)
			outputs = outputs + broadcast_beta

		outputs = K.reshape(outputs, tensor_input_shape)

		return outputs

	def get_config(self):
		config = {'groups': self.groups,'axis': self.axis,'epsilon': self.epsilon,'center': self.center,'scale': self.scale,
			'beta_initializer': initializers.serialize(self.beta_initializer),'gamma_initializer': initializers.serialize(self.gamma_initializer),
			'beta_regularizer': regularizers.serialize(self.beta_regularizer),'gamma_regularizer': regularizers.serialize(self.gamma_regularizer),
			'beta_constraint': constraints.serialize(self.beta_constraint),'gamma_constraint': constraints.serialize(self.gamma_constraint)}
		base_config = super(GroupNormalization, self).get_config()
		return dict(list(base_config.items()) + list(config.items()))

	def compute_output_shape(self, input_shape):
		return input_shape

def UNET(inputs):
	conv1 = Conv1D(32, 6, activation = 'relu', padding = 'same',kernel_regularizer = regularizers.l2(1e-4), bias_regularizer = regularizers.l2(1e-4))(inputs)
	conv1 = GroupNormalization(groups = 4, axis = -1)(conv1)
	pool1 = MaxPooling1D(pool_size=6)(conv1)
	drop1 = Dropout(0.25)(pool1)
	conv2 = Conv1D(64, 3, activation = 'relu', padding = 'same',kernel_regularizer = regularizers.l2(1e-4),  bias_regularizer = regularizers.l2(1e-4))(drop1)
	drop2 = Dropout(0.25)(conv2)
	#pool2 = MaxPooling1D(pool_size=2)(drop2)

	up3 = Conv1D(32, 2, activation = 'relu', padding = 'same', kernel_regularizer = regularizers.l2(1e-4),bias_regularizer = regularizers.l2(1e-4))(drop2)
	drop3 = Dropout(0.25)(up3)
	merge3 = concatenate([drop1,drop3], axis = 2)
	conv4 = Conv1D(32, 3, activation = 'relu', padding = 'same',kernel_regularizer = regularizers.l2(1e-4), bias_regularizer = regularizers.l2(1e-4))(merge3)
	drop4 = Dropout(0.25)(conv4)
	x = Flatten()(drop4)
	x = Dense(16,  kernel_regularizer = regularizers.l2(1e-4),bias_regularizer = regularizers.l2(1e-4))(x)
	x = Activation('relu')(x)
	x = Dropout(0.25)(x)
	#conv3 = Conv1D(16, 1,kernel_regularizer = regularizers.l2(1e-4),bias_regularizer = regularizers.l2(1e-4,activation = 'relu')(conv3)

	return x

def CNN(inputs,filter_size,pool_size):
	conv_initializer = initializers.TruncatedNormal(mean=0., stddev=0.12)
	fc_initializer =   initializers.TruncatedNormal(mean=0., stddev=0.06)
	x = Conv1D(filters= 32, kernel_size= filter_size, padding = 'valid', kernel_regularizer = regularizers.l2(5e-5), bias_regularizer = regularizers.l2(5e-5),kernel_initializer=conv_initializer)(inputs)
	x = GroupNormalization(groups = 4, axis = -1)(x) 
	x = Activation('relu')(x)
	#x = LeakyReLU(alpha=0.1)(x)
	x = MaxPooling1D(pool_size = pool_size)(x)
	x = Dropout(0.25)(x)
	x = Flatten()(x)
	x = Dense(64,  kernel_regularizer = regularizers.l2(5e-4),bias_regularizer = regularizers.l2(5e-4),kernel_initializer=fc_initializer)(x)
	x = Activation('relu')(x)
	#x = LeakyReLU(alpha=0.1)(x)
	x = Dropout(0.25)(x)
	return x

def DCNN(inputs):
	x = Conv1D(filters= 32, kernel_size= 6, padding = 'valid', kernel_regularizer = regularizers.l2(1e-4), bias_regularizer = regularizers.l2(1e-4))(inputs)
	x = GroupNormalization(groups = 4, axis = -1)(x) 
	x = Activation('relu')(x)
	#x = LeakyReLU(alpha=0.1)(x)
	x = MaxPooling1D(pool_size = 6)(x)
	x = Dropout(0.25)(x)
	x = Conv1D(filters= 64, kernel_size= 3, padding = 'valid', kernel_regularizer = regularizers.l2(1e-4),     bias_regularizer = regularizers.l2(1e-4))(x)
	x = Activation('relu')(x)
	x = MaxPooling1D(pool_size = 2)(x)
	x = Dropout(0.25)(x)
	x = Flatten()(x)
	x = Dense(32,  kernel_regularizer = regularizers.l2(1e-4),bias_regularizer = regularizers.l2(1e-4))(x)
	x = Activation('relu')(x)
	#x = LeakyReLU(alpha=0.1)(x)
	x = Dropout(0.25)(x)
	return x

def FC(inputs):
	z = Flatten()(inputs)
	z = Dense(64, kernel_regularizer = regularizers.l2(1e-4),bias_regularizer = regularizers.l2(1e-4))(z)
	#z = Activation('relu')(z)
	z = LeakyReLU(alpha=0.1)(z)
	return z


def PolyA_CNN(length,seq_kernel_size=6,cov_kernel_size=12,pool_size=6):

	input_shape1 = (length,4)
	input_shape2 = (length,1)
	seq_input = Input(shape = input_shape1,name="seq_input")
	cov_input = Input(shape = input_shape2,name="cov_input")
	fc_initializer =   initializers.TruncatedNormal(mean=0., stddev=0.06)

	#
	x = CNN(seq_input,seq_kernel_size,pool_size)
	y = CNN(cov_input,cov_kernel_size,pool_size)
	input_layers = [seq_input,cov_input]
	outLayer= Dense(1, kernel_regularizer = regularizers.l2(5e-4),bias_regularizer = regularizers.l2(5e-4),kernel_initializer=fc_initializer,activation='sigmoid')(concatenate([x,y]))
	model = Model(inputs=input_layers, outputs=outLayer)

	return model
