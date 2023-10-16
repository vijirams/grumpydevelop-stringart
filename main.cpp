/**
 * \file main.cpp
 * \brief implementation of the string art algorithm
 * \author GrumpyDeveloper (Sascha Nitsch)
 * \copyright 2023 Sascha Nitsch
 * Licensed under GPL3 or later license
 * https://contentnation.net/en/grumpydevelop/stringart
 * SPDX-FileCopyrightText: 2023 Sascha Nitsch (@grumpydevelop@contentnation.net) https://contentnation.net/en/grumpydevelop
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */
// compile with
/// g++ -march=native `Magick++-config --cxxflags --cppflags` -Wall -Werror -o main main.cpp -O3 `Magick++-config --ldflags --libs`

#include <math.h>
#include <Magick++.h>
#include <string.h>
#include <sys/time.h>
#include <vector>
#include <unordered_map>

// which basic nail placement algorithms should be used
// #define grid
// #define multicircle
#define circle

/// the actual weight function to calculate how off we are to the target
/// \param value current value
/// \param target desired target
/// \retval distance to target
inline int64_t weightFunction(int16_t value, int16_t target) {
  return (value - target) * (value - target);
}

/// definition of a point for our dwarn line vector
struct Point {
  /// x coordinate
  uint16_t x;
  /// y coordinate
  uint16_t y;
  /// color value
  uint8_t color;
};

/// typedef for a list of points tha make a line
typedef std::vector<Point> td_pointsInLine;
/// the line from src to dst, key = (src << 16) + dst
typedef std::unordered_map<uint32_t, td_pointsInLine> td_linesFromSource;

/// swaps two numbers
/// \param a first number
/// \param b second number
inline void swap(int16_t* a , int16_t* b) {
    int16_t temp = *a;
    *a = *b;
    *b = temp;
}

/// return floating part of number
/// \param x number to process
/// \retval the data behind the dot
inline float fPartOfNumber(float x) {
    if (x > 0) {
      return x - floor(x);
    }
    return x - (floor(x) + 1);
}
/// add given point to vector if col is > 0
/// \param pil pointer to vector
/// \param x x coordinate
/// \param y y coordingte
/// \param col color
inline void pb(td_pointsInLine* pil, uint16_t x, uint16_t y, uint8_t col) {
  if (col > 0) {
    pil->push_back({x, y, col});
  }
}

/// draw line using  Xiaolin Wu’s line algorithm
/// \param x0 source x coordinate
/// \param y0 source y coordinate
/// \param x1 destination x coordinate
/// \param y1 destination y coordinate
/// \param color color to draw
td_pointsInLine drawAALine(int16_t x0 , int16_t y0 , int16_t x1 , int16_t y1, uint8_t color) {
  td_pointsInLine pil;
  bool steep = abs(y1 - y0) > abs(x1 - x0);
  // swap the co-ordinates if slope > 1 or we
  // draw backwards
  if (steep) {
    swap(&x0, &y0);
    swap(&x1, &y1);
  }
  if (x0 > x1) {
    swap(&x0, &x1);
    swap(&y0, &y1);
  }

  // compute the slope
  float dx = x1 - x0;
  float dy = y1 - y0;
  float gradient = dy / dx;
  if (dx == 0.0) {
    gradient = 1;
  }
  int16_t xpxl1 = x0;
  int16_t xpxl2 = x1;
  float intersectY = y0;

  // main loop
  if (steep) {
    int16_t x;
    for (x = xpxl1 ; x <= xpxl2 ; ++x) {
      // pixel coverage is determined by fractional
      // part of y co-ordinate
      pb(&pil, static_cast<uint16_t>(intersectY), static_cast<uint16_t>(x), static_cast<uint8_t>(color * (1 - fPartOfNumber(intersectY))));
      if (intersectY >= 1) {
        pb(&pil, static_cast<uint16_t>(intersectY - 1), static_cast<uint16_t>(x), static_cast<uint8_t>(color * fPartOfNumber(intersectY)));
      }
      intersectY += gradient;
    }
  } else {
    int16_t x;
    for (x = xpxl1 ; x <= xpxl2 ; ++x) {
      // pixel coverage is determined by fractional
      // part of y co-ordinate
      pb(&pil, static_cast<uint16_t>(x), static_cast<uint16_t>(intersectY), static_cast<uint8_t>(color * (1 - fPartOfNumber(intersectY))));
      if (intersectY >= 1) {
        pb(&pil, static_cast<uint16_t>(x), static_cast<uint16_t>(intersectY - 1), static_cast<uint8_t>(color * fPartOfNumber(intersectY)));
      }
      intersectY += gradient;
    }
  }
  return pil;
}

/// \brief main entry point
/// \param argc number of command line arguments
/// \param argv command line arguments
int main(int argc, char* argv[]) {
  if (argc != 8) {
    printf("usage: %s <image name> <resolution x> <resolution y> <number of nails> <max number of iterations> <penalty for duplicate path usage> <lineColor>\n", argv[0]);
    return 1;
  }
  // copy command line data to easier to use variables
  const char* imageName = argv[1];
  uint16_t resolutionX = atoi(argv[2]);
  uint16_t resolutionY = atoi(argv[3]);
  uint16_t numberOfNails = atoi(argv[4]);
  uint16_t maxIter = atoi(argv[5]);
  float duplicateFactor = atof(argv[6]);
  uint8_t lineColor = atoi(argv[7]);

  // our line storage
  td_linesFromSource linesFromSource;

  printf("res: %ix%i nails: %i maxIter: %i duplicatePenalty %.1f color: %i\n", resolutionX, resolutionY, numberOfNails, maxIter, duplicateFactor, lineColor);

  /// nail positions (x << 16) + y
  std::vector<uint32_t> nails;
  // for time measurement
  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  // initialize image magick, load image and resize to target coordinates
  Magick::InitializeMagick(NULL);
  Magick::Image img(imageName);
  if (img.depth() != 8) {
    printf("only 8 bit images supported\n");
    return 1;
  }
  img.sample(Magick::Geometry(resolutionX, resolutionY));

  // fix potential size differences between requested and delivered size
  uint16_t realWidth = img.columns();
  uint16_t realHeight = img.rows();

  // position nails
#ifdef circle
  for (uint16_t i = 0; i < numberOfNails; ++i) {
    float x = sin(2.0 * M_PI * i / numberOfNails) * (realWidth-1) / 2.0 + realWidth / 2.0;
    float y = cos(2.0 * M_PI * i / numberOfNails) * (realHeight-1) / 2.0 + realHeight / 2.0;
    nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
#endif
#ifdef multicircle
  uint16_t count = numberOfNails/1.5;
  for (uint16_t i = 0; i < count; ++i) {
    float x = sin(2.0 * M_PI * i / count) * (realWidth-1) / 2.0 + realWidth / 2.0;
    float y = cos(2.0 * M_PI * i / count) * (realHeight-1) / 2.0 + realHeight / 2.0;
    nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
  uint16_t width = realWidth/1.2 -1;
  uint16_t height = realHeight/1.2 -1;
  count = numberOfNails/1.5;
  for (uint16_t i = 0; i < count; ++i) {
      float x = sin(2.0 * M_PI * i / count) * (width-1) / 2.0 + realWidth / 2.0;
      float y = cos(2.0 * M_PI * i / count) * (height-1) / 2.0 + realHeight / 2.0;
      nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
  width = realWidth/1.5 -1;
  height = realHeight/1.5 -1;
  count = numberOfNails/2;
  for (uint16_t i = 0; i < count; ++i) {
      float x = sin(2.0 * M_PI * i / count) * (width-1) / 2.0 + realWidth / 2.0;
      float y = cos(2.0 * M_PI * i / count) * (height-1) / 2.0 + realHeight / 2.0;
      nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
  width = realWidth/2 -1;
  height = realHeight/2 -1;
  count = numberOfNails/3;
  for (uint16_t i = 0; i < count; ++i) {
      float x = sin(2.0 * M_PI * i / count) * (width-1) / 2.0 + realWidth / 2.0;
      float y = cos(2.0 * M_PI * i / count) * (height-1) / 2.0 + realHeight / 2.0;
      nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
  width = realWidth/3 -1;
  height = realHeight/3 -1;
  count = numberOfNails/4;
  for (uint16_t i = 0; i < count; ++i) {
      float x = sin(2.0 * M_PI * i / count) * (width-1) / 2.0 + realWidth / 2.0;
      float y = cos(2.0 * M_PI * i / count) * (height-1) / 2.0 + realHeight / 2.0;
      nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
  width = realWidth/5 -1;
  height = realHeight/5 -1;
  count = numberOfNails/6;
  for (uint16_t i = 0; i < count; ++i) {
      float x = sin(2.0 * M_PI * i / count) * (width-1) / 2.0 + realWidth / 2.0;
      float y = cos(2.0 * M_PI * i / count) * (height-1) / 2.0 + realHeight / 2.0;
      nails.push_back((static_cast<uint32_t>(floor(x)) << 16) +  static_cast<uint16_t>(floor(y)));
  }
  nails.push_back((static_cast<uint32_t>(floor(realWidth / 2.0)) << 16) +  static_cast<uint16_t>(floor(realHeight / 2.0)));
#endif
#ifdef grid
  uint8_t sq_pins = sqrt(numberOfNails);
  float distX = static_cast<float>(realWidth - 1) / (sq_pins - 1);
  float distY = static_cast<float>(realHeight - 1) / (sq_pins - 1);
  for (uint16_t y = 0; y < sq_pins; ++y) {
    for (uint16_t x = 0; x < sq_pins; ++x) {
      nails.push_back((static_cast<uint32_t>(floor(distX * x)) << 16) +  static_cast<uint16_t>(floor(distY * y)));
    }
  }
#endif
  // number of nails might have been changed above or loaded (in the future)
  numberOfNails = nails.size();
  printf("num %i\n", numberOfNails);

  for (uint16_t src = 0; src < numberOfNails; ++src) {
    for (uint16_t dst = src + 1; dst < numberOfNails; ++dst) {
      td_pointsInLine pointsInLine = drawAALine(nails[src] >> 16, nails[src] & 0xFFFF, nails[dst] >> 16, nails[dst] & 0xFFFF, lineColor);
      linesFromSource.insert(std::make_pair((src << 16) + dst, pointsInLine));
    }
  }

  uint32_t channels = img.channels();
  printf("target image %i x %i x %i\n", realWidth, realHeight, channels);
  MagickCore::Quantum *pixels = img.getPixels(0, 0, realWidth, realHeight);
  uint8_t* targetState = reinterpret_cast<uint8_t*>(malloc(realWidth * realHeight));

  if (channels == 1) {  // monochrome image
    for (uint32_t i = 0; i < realWidth * realHeight; ++i) {
      targetState[i] = pixels[i] >> 8;
    }
  } else if (channels == 2) {  // color + alpha?
    for (uint32_t i = 0; i < realWidth * realHeight; ++i) {
      targetState[i] = pixels[i << 1] >> 8;
    }
  } else {  // RGB or RGBA
    for (uint32_t i = 0; i < realWidth * realHeight; ++i) {
      targetState[i] = ((pixels[i*channels] + pixels[i*channels + 1] + pixels[i*channels + 2]) / 3) >> 8;
    }
  }
#ifdef DEBUGIMG
  FILE* debugfh = fopen("debug.pnm", "wb");
  fprintf(debugfh, "P5\n# debug\n%i %i\n255\n", realWidth, realHeight);
  fwrite(targetState, 1, realWidth * realHeight, debugfh);
  fclose(debugfh);
#endif

  // thread path
  std::vector<uint16_t> path;
  // add start position
  path.push_back(0);
  // a lookup of used paths to count repeats
  uint16_t usedpaths[numberOfNails][numberOfNails];
  bzero(usedpaths, numberOfNails * numberOfNails * 2);

  /// last thread end position
  int16_t lastPosition = 0;
  /// storage for the current state (all previous drawn threads)
  int16_t* currentState = reinterpret_cast<int16_t*>(malloc(realWidth * realHeight * 2));
  /// temp storage to save (current) best version, will be continously updated
  int16_t* bestState = reinterpret_cast<int16_t*>(malloc(realWidth * realHeight * 2));

  // clear states
  uint32_t widthXheight = realWidth * realHeight;
  for (uint32_t i = 0; i < widthXheight; ++i) {
    currentState[i] = 255;
    bestState[i] = 255;
  }
  // current iteration
  uint32_t iter = 0;
  // list of used nails with their counter
  uint8_t usedPins[numberOfNails] = {0};
  // number of continous jump tries if we got stuck
  uint16_t jumps = 0;
  /// total diff from currentState to targetState
  int64_t totalDiff = 0;

  // calculate inital difference
  for (uint32_t i = 0; i < widthXheight; ++i) {
    totalDiff += weightFunction(255, targetState[i]);
  }
  printf("start %li\n", totalDiff);
  while ((iter < maxIter) && jumps*2 < numberOfNails) {
    ++iter;
#ifdef SANITYCHECK
    int64_t sanity = 0;
    for (uint32_t Z = 0; Z < widthXheight; ++Z) {
      int16_t cur = currentState[Z];
      int16_t goal = targetState[Z];
      sanity += weightFunction(cur, goal);
    }
    if (sanity != totalDiff) {
      printf("%i: total: %li, sanity: %li, diff: %li\n", iter, totalDiff, sanity, sanity - totalDiff);
    }
#endif
    /// current best difference
    int64_t realBestDiff = INT64_MAX;
    /// compensated diff includes penality when reusing paths
    int64_t compensatedBestDiff = INT64_MAX;
    int16_t bestTarget = -1;
    // printf("source %i\n", lastPosition); fflush(stdout);
    for (int16_t target = 0; target < numberOfNails; ++target) {
      if (target == lastPosition) continue;
      /// the diff on current lastPosition -> target
      int64_t testDiff = 0;
      uint16_t src = std::min(lastPosition, target);
      uint16_t dst = std::max(lastPosition, target);
      td_linesFromSource::const_iterator lttIter = linesFromSource.find((src << 16) + dst);
      // calculate difference to target
      // for each point
      td_pointsInLine::const_iterator pilIter = lttIter->second.begin();
      while (pilIter != lttIter->second.end()) {
        uint16_t x = (*pilIter).x;
        uint16_t y = (*pilIter).y;
        uint8_t sub = (*pilIter).color;
        uint32_t index = y * realWidth + x;
        int16_t cur = currentState[index];
        int16_t goal = targetState[index];
        // subtract previous error
        testDiff -= weightFunction(cur, goal);
        cur -= sub;
        // add new error
        testDiff += weightFunction(cur, goal);
        ++pilIter;
      }
      float duplicatePenalty = (duplicateFactor != 1) ?
          pow(duplicateFactor, usedpaths[std::min(lastPosition, target)][std::max(lastPosition, target)])
          : 1;
      if ((testDiff / duplicatePenalty) < compensatedBestDiff) {
        // printf("    new best %i - %i(%i,%i) %li d %li\n", lastPosition, target, pos[target]>>16, pos[target]&0xffff, testDiff, totalDiff - testDiff);
        compensatedBestDiff = testDiff * duplicatePenalty;
        realBestDiff = testDiff;
        // a new best
        if (bestTarget != -1) {
          // printf("undo %i:%i\n", lastPosition, bestTarget);
          // undo previous best
          uint16_t tmpSrc = std::min(lastPosition, bestTarget);
          uint16_t tmpDst = std::max(lastPosition, bestTarget);
          td_linesFromSource::const_iterator tmpLfsIter = linesFromSource.find((tmpSrc << 16) + tmpDst);
          td_pointsInLine::const_iterator pilIter = tmpLfsIter->second.begin();
          while (pilIter != tmpLfsIter->second.end()) {
            uint16_t x = (*pilIter).x;
            uint16_t y = (*pilIter).y;
            uint32_t index = y * realWidth + x;
            int16_t cur = currentState[index];
            bestState[index] = cur;
            ++pilIter;
          }
        }
        // apply current best
        td_pointsInLine::const_iterator pilIter = lttIter->second.begin();
        while (pilIter != lttIter->second.end()) {
          uint16_t x = (*pilIter).x;
          uint16_t y = (*pilIter).y;
          int16_t sub = (*pilIter).color;
          uint32_t index = y * realWidth + x;
          int16_t cur = currentState[index];
          cur -= sub;
          bestState[index] = cur;
          ++pilIter;
        }
        bestTarget = target;
      }
    }
    if (realBestDiff >= 0) {
      // we got worse, jump to random place to continue
      // printf("j %3i %3i -> %3i(%4i, %4i) bestDiff %8li (%li) iter %i path %i\n", jumps, lastPosition, bestTarget, pos[bestTarget] >> 16, pos[bestTarget] & 0xFFFF, realBestDiff, totalDiff - realBestDiff, iter, usedpaths[std::min(lastPosition, bestTarget)][std::max(lastPosition, bestTarget)]);
      if (jumps) {  // undo last jump, was not working anyway
        --usedPins[lastPosition];
        path.pop_back();
      }
      // select next target randomly (kind of, intentially producing the same numbers)
      bestTarget = random() % numberOfNails;
      path.push_back(bestTarget);
      lastPosition = bestTarget;
      ++usedPins[bestTarget];
      ++jumps;
      --iter;
      // fix bestState, easier to copy over than manually reversing. should only happen a few times anyway
      memcpy(bestState, currentState, widthXheight * 2);
      continue;
    }
    /// reset jump counter (faster to always set)
    jumps = 0;
    if (bestTarget < 0) {
      printf("no best\n");
      break;
    }
    // update path map
    ++usedpaths[std::min(lastPosition, bestTarget)][std::max(lastPosition, bestTarget)];
    // add new stop
    path.push_back(bestTarget);
    // update used pins
    ++usedPins[bestTarget];
    // progress report
    if (iter % 100 == 0) {
      printf("best %4i -> %4i(%4i, %4i) diff %9li (%12li) iter %5i path %i\n", lastPosition, bestTarget, nails[bestTarget] >> 16, nails[bestTarget] & 0xFFFF, compensatedBestDiff, totalDiff + realBestDiff, iter, usedpaths[std::min(lastPosition, bestTarget)][std::max(lastPosition, bestTarget)]);
    }
    // set new start position
    lastPosition = bestTarget;
    // update diff
    totalDiff += realBestDiff;
    // update current state from best map
    memcpy(currentState, bestState, widthXheight * 2);
  }

  printf("size %li\n", path.size());
  // we are done, create output svg
  FILE* fh = fopen("map.svg", "wb");
  fprintf(fh, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 %i %i\">\n<rect width=\"%i\" height=\"%i\" fill=\"#ffffff\" />\n<g style=\"fill:none;stroke:#000000;stroke-opacity:%.2f;stroke-width:1\">\n", realWidth, realHeight, realWidth, realHeight, lineColor / 255.0);
  uint32_t counter = 0;
  for (uint16_t i : path) {
    if ((counter & 255) == 0) {
      fprintf(fh, "<path d=\"M%i %i", nails[i] >> 16, nails[i] & 0xffff);
    } else {
      fprintf(fh, "L%i %i", nails[i] >> 16, nails[i] & 0xffff);
    }
    if ((counter & 255) == 255) {
      fprintf(fh, "\" />\n");
    }
    ++counter;
  }
  if ((counter & 255) != 0) {
    fprintf(fh, "\" />\n");
  }
gettimeofday(&tv2, NULL);
float timeNeeded = (tv2.tv_sec - tv1.tv_sec) + (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
  fprintf(fh, "</g>\n");
  fprintf(fh, "<text x=\"30\" y=\"10\" style=\"font-weight:bold;font-size:60px;font-family:'DejaVu Serif'\"><tspan x=\"30\" y=\"70\">%s %s %i %i</tspan><tspan x=\"70\" y=\"150\"> %i %i %.2f %i</tspan><tspan x=\"30\" y=\"%i\">%i nails, %li paths, %.1f sec</tspan></text>", argv[0], imageName, resolutionX, resolutionY, atoi(argv[4]), maxIter, duplicateFactor, lineColor, realHeight - 30, numberOfNails, path.size(), timeNeeded);
  fprintf(fh, "</svg>");
  fclose(fh);
  // cleanup
  free(targetState);
  free(bestState);
  free(currentState);
  path.clear();
  linesFromSource.clear();
  Magick::TerminateMagick();
}
