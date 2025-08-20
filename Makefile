CC = gcc
CFLAGS = -g -Wall
LDFLAGS = -lm
SRCS = dpsk_fec.c reedsolomon.c fft-complex.c
OBJS = $(SRCS:.c=.o)
TARGET = dpsk_fec

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
