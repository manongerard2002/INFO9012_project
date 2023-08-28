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

#include <mpi.h>

#include "colors.h"

float futureTime;
int width;
int height;

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

int main(int argc, char **argv)
{
  int rank, size, i;
  long double sum = 0, partialSum = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int futureTime;
  double scale, cx, cy, angle;  

  // became local
  scale = 1./256;
  cx = -0.74364389984420875557,
  cy = 0.13182603113504143266;
  angle = 0;
  float fps = 1.;
  unsigned keys = 0;
  bool rotate = false;
  timer_reset();
  
  if (argc < 4)
    goto usage;

  // became global
  width = atoi(argv[1]);
  height = atoi(argv[2]);
  int nthreads = atoi(argv[3]); //not needed
  int sizePdata = width * height;

  unsigned *Pdata;
  Pdata = malloc(sizePdata * 4);
  if (!Pdata)
    handle_error("malloc error");
if (rank == 0) { // the master/consumer
  // 1 consumer:
  // - displaying images,
  // - dealing with user interactions and
  // - require images to be computed
  if (width < 50 || height < 50 || nthreads < 1 || nthreads > 16)
    goto usage;

  Display *display;
  if (!(display = XOpenDisplay(NULL)))
    handle_error("Unable to open display");
  int screen = XDefaultScreen(display);
  if (ScreenOfDisplay(display, screen)->width < width ||
      ScreenOfDisplay(display, screen)->height < height)
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

 skip_usage:  ; //This is an empty statement.
  unsigned long black = BlackPixel(display, screen);
  unsigned long white = WhitePixel(display, screen);

  Window win = XCreateSimpleWindow(display, DefaultRootWindow(display), 0, 0,
                                   width, height, 5, white, black);

  Visual *visual = DefaultVisual(display, screen);

  XImage *img = XCreateImage(display, visual,
                             DefaultDepth(display, screen), ZPixmap, 0,
                             (char *)Pdata, width, height, 32, 0);
  if (!img)
    handle_error("error creating image");
  XMapWindow(display, win);
  GC gc = XCreateGC(display, win, 0, NULL);
  XSelectInput(display, win,
               ExposureMask | KeyPressMask | KeyReleaseMask | 
               ButtonPressMask | StructureNotifyMask);

  printf("FPS: %5.1f", fps);

  // setting the first images (in number equal to the number of computing threads)
  // to compute with a “to compute” status at an interval of 1 second
  // with the appropriate parameters (scale, x, y, angle)
  // future time décrit le temps entre 2 frames (actuel et celle d'avant)
  futureTime = 1000;
  int i;
  for (i = 1; i < size; i++)
  {
    MPI_Send(&futureTime, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    MPI_Send(&scale, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    MPI_Send(&cx, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    MPI_Send(&cy, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    MPI_Send(&angle, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
  }
  long ms = 0;
  float fps_target = 1; // fps voulu initialisé a 1
  timer_reset();
  int time;

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
      i = (i + 1) % size;
      if (i == 0)
        i = 1;
      //dequeue the image
      MPI_Recv(Pdata, sizePdata, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&time, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      XPutImage(display, win, gc, img, 0, 0, 0, 0, width, height);
      XFlush(display);
      if(timer_get() <= time)
      {
        if (fps_target*1.1 >= 30)
        {
          fps_target = 30; // fps is capped to 30
        }
        else
          fps_target = fps_target * 1.1;
        while(timer_get() <= time)
        {
          usleep((time-timer_get())*250);
        }
      }
      else{
        fps_target = fps_target / 1.1;
      }

      futureTime = (long) (1000.0/fps_target);
      //enqueue the task
      MPI_Send(&futureTime, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&scale, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      MPI_Send(&cx, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      MPI_Send(&cy, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      MPI_Send(&angle, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

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
  XFreeGC(display, gc);
  XCloseDisplay(display);
  }
  else { // producer
  // send msg to master
    while (true)
    {
      // look for work to do : wait receiving parameters
      MPI_Recv(&futureTime, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&scale, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&cx, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&cy, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&angle, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // do the work : update Pdata
      mandelbrot(Pdata, scale, cx, cy, angle);

      // Send the Pdata and future time
      MPI_Send(Pdata, sizePdata, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&futureTime, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }
  free(Pdata);
  MPI_Finalize();

  return 0;
}
