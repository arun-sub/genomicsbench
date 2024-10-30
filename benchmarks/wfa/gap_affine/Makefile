###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT:=..
FOLDER_BUILD_PATH:=$(FOLDER_ROOT)/$(FOLDER_BUILD)

###############################################################################
# Modules
###############################################################################
MODULES=affine_wavefront \
        affine_wavefront_align \
        affine_wavefront_backtrace \
        affine_wavefront_extend \
        affine_wavefront_penalties \
        affine_wavefront_reduction \
        affine_wavefront_utils \
        edit_cigar

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD_PATH)/, $(SRCS:.c=.o))

CC_XFLAGS=$(ARCH_FLAGS)

###############################################################################
# Rules
###############################################################################
        
all: $(OBJS)

$(FOLDER_BUILD_PATH)/affine_wavefront.o : affine_wavefront.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD_PATH)/affine_wavefront_align.o : affine_wavefront_align.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD_PATH)/affine_wavefront_extend.o : affine_wavefront_extend.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
# General building rule
$(FOLDER_BUILD_PATH)/%.o : %.c
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	

