/*
 * gcc -Wall -W -g3 -o mandelb mandelb.c -lX11 -lm
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <sys/time.h>

#include <pthread.h>
#include <error.h>

#include "colors.h"

float futureTime;
pthread_mutex_t lock_buffer_0;
pthread_mutex_t lock_buffer_1;
pthread_mutex_t lock_buffer_2;
pthread_mutex_t lock_buffer_3;
pthread_mutex_t lock_buffer_4;
pthread_mutex_t lock_buffer_5;
pthread_mutex_t lock_buffer_6;
pthread_mutex_t lock_buffer_7;
pthread_mutex_t lock_buffer_8;
pthread_mutex_t lock_buffer_9;
pthread_mutex_t lock_buffer_10;
pthread_mutex_t lock_buffer_11;
pthread_mutex_t lock_buffer_12;
pthread_mutex_t lock_buffer_13;
pthread_mutex_t lock_buffer_14;
pthread_mutex_t lock_buffer_15;
int width;
int height;

pthread_mutex_t chooseMutex(unsigned idx)
{
  switch (idx)
  {
    case 0:
      return lock_buffer_0;
    case 1:
      return lock_buffer_1;
    case 2:
      return lock_buffer_2;
    case 3:
      return lock_buffer_3;
    case 4:
      return lock_buffer_4;
    case 5:
      return lock_buffer_5;
    case 6:
      return lock_buffer_6;
    case 7:
      return lock_buffer_7;
    case 8:
      return lock_buffer_8;
    case 9:
      return lock_buffer_9;
    case 10:
      return lock_buffer_10;
    case 11:
      return lock_buffer_11;
    case 12:
      return lock_buffer_12;
    case 13:
      return lock_buffer_13;
    case 14:
      return lock_buffer_14;
    default:
      return lock_buffer_15;
  }
}

void LockBufferElement(unsigned idx)
{
  switch (idx)
  {
    case 0:
      pthread_mutex_lock(&lock_buffer_0);
      break;
    case 1:
      pthread_mutex_lock(&lock_buffer_1);
      break;
    case 2:
      pthread_mutex_lock(&lock_buffer_2);
      break;
    case 3:
      pthread_mutex_lock(&lock_buffer_3);
      break;
    case 4:
      pthread_mutex_lock(&lock_buffer_4);
      break;
    case 5:
      pthread_mutex_lock(&lock_buffer_5);
      break;
    case 6:
      pthread_mutex_lock(&lock_buffer_6);
      break;
    case 7:
      pthread_mutex_lock(&lock_buffer_7);
      break;
    case 8:
      pthread_mutex_lock(&lock_buffer_8);
      break;
    case 9:
      pthread_mutex_lock(&lock_buffer_9);
      break;
    case 10:
      pthread_mutex_lock(&lock_buffer_10);
      break;
    case 11:
      pthread_mutex_lock(&lock_buffer_11);
      break;
    case 12:
      pthread_mutex_lock(&lock_buffer_12);
      break;
    case 13:
      pthread_mutex_lock(&lock_buffer_13);
      break;
    case 14:
      pthread_mutex_lock(&lock_buffer_14);
      break;
    default:
      pthread_mutex_lock(&lock_buffer_15);
      break;
  }
}

void UnlockBufferElement(unsigned idx)
{
  switch (idx)
  {
    case 0:
      pthread_mutex_unlock(&lock_buffer_0);
      break;
    case 1:
      pthread_mutex_unlock(&lock_buffer_1);
      break;
    case 2:
      pthread_mutex_unlock(&lock_buffer_2);
      break;
    case 3:
      pthread_mutex_unlock(&lock_buffer_3);
      break;
    case 4:
      pthread_mutex_unlock(&lock_buffer_4);
      break;
    case 5:
      pthread_mutex_unlock(&lock_buffer_5);
      break;
    case 6:
      pthread_mutex_unlock(&lock_buffer_6);
      break;
    case 7:
      pthread_mutex_unlock(&lock_buffer_7);
      break;
    case 8:
      pthread_mutex_unlock(&lock_buffer_8);
      break;
    case 9:
      pthread_mutex_unlock(&lock_buffer_9);
      break;
    case 10:
      pthread_mutex_unlock(&lock_buffer_10);
      break;
    case 11:
      pthread_mutex_unlock(&lock_buffer_11);
      break;
    case 12:
      pthread_mutex_unlock(&lock_buffer_12);
      break;
    case 13:
      pthread_mutex_unlock(&lock_buffer_13);
      break;
    case 14:
      pthread_mutex_unlock(&lock_buffer_14);
      break;
    default:
      pthread_mutex_unlock(&lock_buffer_15);
      break;
  }
}

/* Maximum amount of iterations */
#define MAX_ITER 2000

#define BAILOUT 128
#define Q1LOG2 1.44269504088896340735992468100189213742664595415299
#define LOG2 0.69314718055994530941723212145817656807550013436026
#define LOGLOGBAILOUT log(log(BAILOUT))

#define KEY_UP    (1 << 0)
#define KEY_DOWN  (1 << 1)
#define KEY_LEFT  (1 << 2)
#define KEY_RIGHT (1 << 3)
#define KEY_SPACE (1 << 4)
#define KEY_MINUS (1 << 5)
#define KEY_B     (1 << 6)
#define KEY_R     (1 << 7)
/**
   \brief print an error message and exit
   \param message the message to print */
static void
handle_error(char * message)
{
  fprintf(stderr, "%s\n", message);
  exit(EXIT_FAILURE);
}

/** \brief timer facility */ 
static struct timeval timer_start;

/** \brief reset the timer */
static void
timer_reset(void)
{
  gettimeofday(&timer_start, NULL);
}

/** \brief returns the value of the timer in ms */
static long
timer_get(void)
{
  struct timeval now;
  gettimeofday(&now, NULL);
  return (now.tv_sec - timer_start.tv_sec) * 1000 +
    (now.tv_usec - timer_start.tv_usec) / 1000;
}

#if 0

/* This might be useful to get the absolute time: use time_start at the
   beginning of your program, and time_get will provide the number of ms
   since the start */

static struct timeval time_zero;

static void
time_start(void)
{
  gettimeofday(&time_zero, NULL);
}

/* In principle, one must be careful of overflow, but we are
   certainly safe here */

static long
time_get(void)
{
  struct timeval now;
  gettimeofday(&now, NULL);
  return (now.tv_sec - time_zero.tv_sec) * 1000 +
    (now.tv_usec - time_zero.tv_usec) / 1000;
}

#endif

/**
   \brief convert a double value into a byte
   \param d double the value in [0,1[ */
#define to_byte(d) (((int) (d * 256)) & 255)

#if 0
/**
   \brief advanced color function for mandelbrot
   \param iter number of iteration
   \param z2 square of modulus at escape */
static inline unsigned
rgb(unsigned iter, double z2)
{
  double v = log(iter + 1.5 - log(log(z2) / 2.) / log(2.)) / 3.4;
  if (v < 1.0)
    return RGB(to_byte(pow(v, 4)), to_byte(pow(v, 2.5)), to_byte(v));
  v = 2.0 - v;
  if (v < 0.) v = 0.;
  return RGB(to_byte(v), to_byte(pow(v, 1.5)), to_byte(pow(v, 3)));
}
#endif

/**
   \brief advanced color function for mandelbrot
   \param iter number of iteration
   \param z2 square of modulus at escape */
static inline unsigned
rgb2(unsigned iter, double z2)
{
  long double c = iter - 1.28 + (LOGLOGBAILOUT - log(log(z2) / 2.)) * Q1LOG2;
  unsigned idx = fmod((log(c/64+1)/LOG2+0.45), 1) * GRADIENTLENGTH + 0.5;
  return colors[idx];
}

/**
   \brief computes the mandelbrot figure
   \param Pdata the buffer containing the image
   \param scale the scale of the figure
   \param cx the x coordinate of the center point
   \param cy the y coordinate of the center point
   \param angle the angle of the figure
   \note other (global) parameters are width (the width of the image) and height (the height of the image)
 */
static void
mandelbrot(unsigned *Pdata, double scale, double cx, double cy, double angle)
{
  for (int i = 0; i < height; i++) {
    double delta_y = (i - height / 2) * scale;
    for (int j = 0; j  < width; j++, Pdata++) {
      double delta_x = (j - width / 2) * scale;
      double y = delta_y * cos(angle) - delta_x * sin(angle) + cy;
      double x = delta_x * cos(angle) + delta_y * sin(angle) + cx;
      double zx = 0, zy = 0, zx2 = 0, zy2 = 0;
      unsigned iter = 0;
      for (; iter < MAX_ITER && zx2 + zy2 < BAILOUT; iter++)
        {
          zy = 2 * zx * zy + y;
          zx = zx2 - zy2 + x;
          zx2 = zx * zx;
          zy2 = zy * zy;
        }
      *Pdata = (iter >= MAX_ITER) ? 0 : rgb2(iter, zx2+zy2);
    }
  }
}

typedef enum
{
  compute,
  computing,
  computed,
  unassigned,
} StatusType;

typedef struct Element_t
{
  int time;
  StatusType status;
  unsigned *Pdata; // image buffer
  XImage *image; // XImage
  // parameters :
  double scale;
  double cx;
  double cy;
  double angle;
} Element;

// A structure representing a circular buffer.
typedef struct
{
  unsigned i_P;      // write "compute" pointer
  unsigned i_U;      // write "unassigned" pointer
  unsigned i_C;      // read "computed" pointer
  unsigned capacity; // capacity of circular buffer //in reality - 1
  Element *buffer;   // table of Element ?
} circularBuffer;

/*Prototypes*/
circularBuffer *circularBufferCreate(int capacity, Display *display, Visual *visual, int screen);
bool circularBufferIsEmpty(circularBuffer *cbuf);
bool circularBufferIsFull(circularBuffer *cbuf);
bool circularBufferEnqueueToCompute(circularBuffer *cbuf, int futureTime, double scale, double cx, double cy, double angle);
bool circularBufferChooseNextToCompute(circularBuffer *cbuf, unsigned *index);
bool circularBufferUpdateComputed(circularBuffer *cbuf, unsigned idx);
bool circularBufferDequeue(circularBuffer *cbuf, unsigned *index);
void circularBufferFree(circularBuffer *cbuf);
void bufferElementFree(Element *bufferElement);
int bufferGetFutureTime(Element *bufferElement);
XImage *bufferGetImage(Element *bufferElement);
 
/* Create an empty circularBuffer, with all tasks with a status unassigned.
 * The associated Pdata and images are also allocated
 *
 * PARAMETERS:
 * - `capacity` : the size of the circularBuffer
 * - `width` : the width of the window
 * - `height` : the height of the window
 *
 * RETURNS: the new circularBuffer.
 */
circularBuffer *circularBufferCreate(int capacity, Display *display, Visual *visual, int screen)
{
  circularBuffer *cbuf = malloc(sizeof(circularBuffer));
  if (!cbuf)
    handle_error("malloc error");

  cbuf->i_P = 0;
  cbuf->i_U = 0;
  cbuf->i_C = 0;
  cbuf->capacity = capacity;
  cbuf->buffer = malloc(capacity * sizeof(Element));
  if (!cbuf->buffer)
  {
    free(cbuf);
    handle_error("malloc error");
  }
  for (int i = 0; i < capacity; i++)
  {
    cbuf->buffer[i].status = unassigned;
    cbuf->buffer[i].Pdata = malloc(width * height * 4);
    if (!cbuf->buffer[i].Pdata)
    {
      for (int j = 0; j < i; j++)
        free(cbuf->buffer[j].Pdata);
      free(cbuf->buffer);
      free(cbuf);
      handle_error("malloc error");
    }
    cbuf->buffer[i].image = XCreateImage(display, visual, DefaultDepth(display, screen), ZPixmap, 0, (char *)cbuf->buffer[i].Pdata, width, height, 32, 0);
    if (!cbuf->buffer[i].image)
    {
      for (int j = 0; j <= i; j++)
        free(cbuf->buffer[j].Pdata);
      free(cbuf->buffer);
      free(cbuf);
      handle_error("error creating image");
    }
  }

  return cbuf;
}

/* Check whether or not a circularBuffer is empty.
 * Empty : i_U == i_C
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 *
 * RETURNS: true if the circularBuffer is empty, false otherwise.
 */
bool circularBufferIsEmpty(circularBuffer *cbuf)
{
  return (cbuf->i_C == cbuf->i_U);
}

/* Check whether or not a circularBuffer is full.
 * Full : (i_U + 1) % capacity == i_C
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 *
 * RETURNS: true if the circularBuffer is full, false otherwise.
 */
bool circularBufferIsFull(circularBuffer *cbuf)
{
  return (((cbuf->i_U + 1) % cbuf->capacity) == cbuf->i_C);
}

/* Enqueue a new element "to compute" on a circularBuffer, if not full.
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 * - `futureTime`: the time to push on the circularBuffer.
 * - `scale`: the scale of the figure
 * - `cx`: the x coordinate of the center point
 * - `cy`: the y coordinate of the center point
 * - `angle`: the angle of the figure
 *
 * RETURNS: true if not full/succeeded, false otherwise/full/error.
 */
bool circularBufferEnqueueToCompute(circularBuffer *cbuf, int futureTime, double scale, double cx, double cy, double angle)
{
  unsigned idx = cbuf->i_U;
  LockBufferElement(idx);
  bool isFull = circularBufferIsFull(cbuf);
  if (!isFull)
  {
    cbuf->buffer[idx].time = futureTime;
    cbuf->buffer[idx].status = compute;
    cbuf->buffer[idx].scale = scale;
    cbuf->buffer[idx].cx = cx;
    cbuf->buffer[idx].cy = cy;
    cbuf->buffer[idx].angle = angle;
    cbuf->i_U = (idx + 1) % cbuf->capacity;
  }
  UnlockBufferElement(idx);

  return (!isFull);
}

/* Choose the data that is next to be computed,
 * change its status to computing and return it's element and index in argument.
 * If sucessful, keeps the lock on the element
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 * - `data` : a pointer where to return the data
 * - `index` : a pointer where to return the index
 *
 * RETURNS: true if not empty/succeeded, false otherwise/none to compute.
 */
bool circularBufferChooseNextToCompute(circularBuffer *cbuf, unsigned *index)
{
  unsigned idx = cbuf->i_P;
  LockBufferElement(idx);
  bool isCompute = (cbuf->buffer[idx].status == compute);
  if (isCompute)
  {
    cbuf->buffer[idx].status = computing;
    *index = idx;
    cbuf->i_P = (idx + 1) % cbuf->capacity;
  }
  else
    UnlockBufferElement(idx);

  return isCompute;
}

/* Update an element "computed" on a circularBuffer, if not full.
 * + Unlock the element
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 * - `idx`: the index in the buffer
 * - `newPdata`: the Pdata calculated
 *
 * RETURNS: true if succeeded, false otherwise/not previously computing.
 */
bool circularBufferUpdateComputed(circularBuffer *cbuf, unsigned idx)
{
  bool isComputing = (cbuf->buffer[idx].status == computing);
  if (isComputing)
  {
    cbuf->buffer[idx].status = computed;
  }
  UnlockBufferElement(idx);

  return isComputing;
}

/* Dequeue the data that is computed circularBuffer and return it.
 * If dequeued, keeps the lock on the element
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 * - `data` : a pointer where to return the data
 * - `index` : a pointer where to return the index
 *
 * RETURNS: true if not succeeded, false otherwise/none computed.
 */
bool circularBufferDequeue(circularBuffer *cbuf, unsigned *index)
{
  unsigned idx = cbuf->i_C;
  LockBufferElement(idx);
  bool isComputed = (cbuf->buffer[idx].status == computed);
  if (isComputed)
  {
    cbuf->buffer[idx].status = unassigned;
    cbuf->i_C = (idx + 1) % cbuf->capacity;
    *index = idx;
  }
  else
    UnlockBufferElement(idx);

  return isComputed;
}

/* Deallocate the memory occupied by a circularBuffer.
 *
 * PARAMETERS:
 * - `cbuf`: the circularBuffer
 */
void circularBufferFree(circularBuffer *cbuf)
{
  for (int i = 0; i < cbuf->capacity; i++)
  {
    bufferElementFree(&(cbuf->buffer[i]));
  }
  free(cbuf->buffer);
  free(cbuf);
}

/* Deallocate the memory occupied by an element of the buffer.
 *
 * PARAMETERS:
 * - `bufferElement`: the element of the buffer
 */
void bufferElementFree(Element *bufferElement)
{
  XDestroyImage(bufferElement->image);
}

circularBuffer *cbuf;

// mythread of the producer
void *producer(void *arg)
{
  while (true)
  {
    unsigned idx;

    // look for work to do (to compute) and change the status (computing)
    while (!circularBufferChooseNextToCompute(cbuf,  &idx))
    {
      usleep(16670);
    }

    // do the work : update Pdata
    mandelbrot(cbuf->buffer[idx].Pdata, cbuf->buffer[idx].scale, cbuf->buffer[idx].cx, cbuf->buffer[idx].cy, cbuf->buffer[idx].angle);

    // change the status again (computed)
    circularBufferUpdateComputed(cbuf, idx);
  }
  return NULL;
}

int main(int argc, char **argv)
{
  // became local
  double scale = 1./256;
  double cx = -0.74364389984420875557,
         cy = 0.13182603113504143266;
  double angle = 0;
  float fps = 1.;
  unsigned keys = 0;
  bool rotate = false;
  unsigned *Pdata;
  timer_reset();
  
  if (argc < 4)
    goto usage;

  // became global
  width = atoi(argv[1]);
  height = atoi(argv[2]);
  int nthreads = atoi(argv[3]);

  if (width < 50 || height < 50 || nthreads < 1 || nthreads > 16)
    goto usage;

  Display *display;
  if (!(display = XOpenDisplay(NULL)))
    handle_error("Unable to open display");
  int screen = XDefaultScreen(display);
  if (ScreenOfDisplay(display, screen)->width < width ||
      ScreenOfDisplay(display, screen)->height < height || nthreads < 2)
    {
      XCloseDisplay(display);
      goto usage;
    }
  goto skip_usage;
 usage:
  fprintf(stderr, "usage: %s <width> <height> <nthreads> \n", argv[0]);
  fprintf(stderr, "  <width> and <height> must be greater than or equal to 50\n");
  fprintf(stderr, "  <width> and <height> must be smaller than or equal to the screen size\n");
  fprintf(stderr, "  <nthreads> be an integer in {2,... 16}\n");
  return 1;

 skip_usage:
  Pdata = malloc(width * height * 4);
  if (!Pdata)
    handle_error("malloc error");

  unsigned long black = BlackPixel(display, screen);
  unsigned long white = WhitePixel(display, screen);

  Window win = XCreateSimpleWindow(display, DefaultRootWindow(display), 0, 0,
                                   width, height, 5, white, black);

  Visual *visual = DefaultVisual(display, screen);

  XMapWindow(display, win);
  GC gc = XCreateGC(display, win, 0, NULL);
  XSelectInput(display, win,
               ExposureMask | KeyPressMask | KeyReleaseMask | 
               ButtonPressMask | StructureNotifyMask);

  printf("FPS: %5.1f", fps);

  // create buffer element locks
  for (unsigned i = 0; i < 16; i++)
  {
    pthread_mutex_t lock_buffer = chooseMutex(i);
    if (pthread_mutex_init(&lock_buffer, NULL))
      error(-1, 0, "Error initializing mutex");
  }

  // create the circular queue
  // faudra réflechir a sa capacity par essais erreurs:
  // 16 : pour avoir au moins la place de push les 15 trucs initiaux
  cbuf = circularBufferCreate(16, display, visual, screen);

  // create a thread pool
  // n-1 producers: the proper computation of images
  unsigned i;
  pthread_t thread[nthreads-1];
  for (i = 0; i < nthreads - 1; i++)
  {
    if (pthread_create(&thread[i], NULL, producer, NULL))
      error(-1, 0, "Error creating thread");
  }

  // 1 consumer:
  // - displaying images,
  // - dealing with user interactions and
  // - require images to be computed

  // setting the first images (in number equal to the number of computing threads)
  // to compute with a “to compute” status at an interval of 1 second
  // with the appropriate parameters (scale, x, y, angle)
  // future time décrit le temps entre 2 frames (actuel et celle d'avant)
  futureTime = 1000;
  for (int i = 0; i < nthreads - 1; i++)
  {
    while(!circularBufferEnqueueToCompute(cbuf, futureTime, scale, cx, cy, angle));
  }
  fprintf(stderr, "\n");
  long ms = 0;
  float fps_target = 1; // fps voulu initialisé a 1
  timer_reset();
  unsigned idx;

  while (1)
    {
      XEvent event;
      while (XPending(display))
        {
          XNextEvent(display, &event);
          switch (event.type) {
          case Expose:
            break;            
          case KeyPress:
            switch (XLookupKeysym(&event.xkey, 0)) {
            case XK_Up :    keys |= KEY_UP;    break;
            case XK_Down :  keys |= KEY_DOWN;  break;
            case XK_Left :  keys |= KEY_LEFT;  break;
            case XK_Right : keys |= KEY_RIGHT; break;
            case XK_space : keys |= KEY_SPACE; break;
            case XK_minus : keys |= KEY_MINUS; break;
            case XK_b :     keys |= KEY_B; break;
            case XK_r :
              if (!(keys & KEY_R)) rotate = !rotate;
              keys |= KEY_R; break;
            case XK_q : goto exiting;
            }
            break;
          case KeyRelease:
            switch (XLookupKeysym(&event.xkey, 0)) {
            case XK_Up :    keys &= ~KEY_UP;    break;
            case XK_Down :  keys &= ~KEY_DOWN;  break;
            case XK_Left :  keys &= ~KEY_LEFT;  break;
            case XK_Right : keys &= ~KEY_RIGHT; break;
            case XK_space : keys &= ~KEY_SPACE; break;
            case XK_minus : keys &= ~KEY_MINUS; break;
            case XK_b :     keys &= ~KEY_B; break;
            case XK_r :     keys &= ~KEY_R; break;
            }
            break;
          case ClientMessage:
            goto exiting;
          }
        }
      if (keys & KEY_UP)
        {
          cy -= 1 * cos(angle) * scale;
          cx -= 1 * sin(angle) * scale;
        }
      if (keys & KEY_DOWN)
        {
          cy += 1 * cos(angle) * scale;
          cx += 1 * sin(angle) * scale;
        }
      if (keys & KEY_LEFT)
        {
          cx -= 1 * cos(angle) * scale;
          cy += 1 * sin(angle) * scale;
        }
      if (keys & KEY_RIGHT)
        {
          cx += 1 * cos(angle) * scale;
          cy -= 1 * sin(angle) * scale;
        }
      if (keys & KEY_SPACE)
        scale /= 1.1;
      if ((keys & KEY_MINUS) || (keys & KEY_B))
        scale *= 1.1;
      if (rotate)
        {
          angle += 2 * 3.1415 / fps / 30;
          if (angle > 2 * 3.1415) angle -= 2 * 3.1415;
        }
      // busy wait for available images (computed status)
      while (!circularBufferDequeue(cbuf, &idx))
      {
        usleep(16670);
      }

      XPutImage(display, win, gc, cbuf->buffer[idx].image, 0, 0, 0, 0, width, height);
      XFlush(display);
      if(timer_get() <= cbuf->buffer[idx].time)
      {
        if (fps_target*1.1 >= 30)
        {
          fps_target = 30; // fps is capped to 30
        }
        else
          fps_target = fps_target * 1.1;
        while(timer_get() <= cbuf->buffer[idx].time)
        {
          usleep((cbuf->buffer[idx].time-timer_get())*250);
        }
      }
      else{
        fps_target = fps_target / 1.1;
      }
      
      UnlockBufferElement(idx);
      futureTime = (long) (1000.0/fps_target);
      while(circularBufferEnqueueToCompute(cbuf, futureTime, scale, cx, cy, angle));

      static unsigned n = 0;
      n++;
      ms += timer_get();
      if (n > 10 || (n && ms > 2000))
      {
        fps = n * 1000.0 / ms;
        printf("\rFPS: %4.1f", fps);
        fflush(stdout);
        n = 0;
        ms = 0;
      }
      timer_reset();
    }
 exiting:
  // stop thread before ending the program
  for (i = 0; i < nthreads - 1; i++)
    if (pthread_cancel(thread[i]))
      error(-1, 0, "Error stopping thread");

  // destroy buffer element locks
  for (unsigned i = 0; i < nthreads-1; i++)
  {
    pthread_mutex_t lock_buffer = chooseMutex(i);
    if (pthread_mutex_destroy(&lock_buffer))
      error(-1, 0, "Error destroying mutex");
  }
  circularBufferFree(cbuf);
  XFreeGC(display, gc);
  XCloseDisplay(display);
  return 0;
}
