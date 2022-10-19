/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019 Michael Carley
 *
 * WBFMM is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  WBFMM is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with WBFMM.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @defgroup logging Logging functions
 * 
 * @brief  Logging functions for use with GLIB logging facilities.
 * 
 * The logging functions provide control over the messages logged by
 * WBFMM, to ease debugging and to give access to messages from
 * specific processors during parallel execution.
 *
 * @{
 */

#include <stdio.h>

#include <glib.h>

#include "wbfmm.h"

#define WBFMM_LOGGING_DATA_WIDTH     8
#define WBFMM_LOGGING_DATA_FID       0
#define WBFMM_LOGGING_DATA_PREFIX    1
#define WBFMM_LOGGING_DATA_LEVEL     2
#define WBFMM_LOGGING_DATA_EXIT_FUNC 3
#define WBFMM_LOGGING_DATA_TIMER     4

void wbfmm_logging_func(const gchar *log_domain,
		       GLogLevelFlags log_level,
		       const gchar *message,
		       gpointer data[]) ;

const gchar *wbfmm_logging_string(GLogLevelFlags level) ;

/** 
 * Return a string describing the log level of the message.
 * 
 * @param level log level.
 * 
 * @return string describing level.
 */

const gchar *wbfmm_logging_string(GLogLevelFlags level)

{
  const gchar *strings[] = {"RECURSION", 
			    "FATAL",
			    "ERROR",
			    "CRITICAL",
			    "WARNING",
			    "MESSAGE",
			    "INFO",
			    "DEBUG"} ;

  if ( G_LOG_LEVEL_ERROR & level) return strings[2] ; 
  if ( G_LOG_LEVEL_CRITICAL & level) return strings[3] ; 
  if ( G_LOG_LEVEL_WARNING & level) return strings[4] ; 
  if ( G_LOG_LEVEL_MESSAGE & level) return strings[5] ; 
  if ( G_LOG_LEVEL_INFO & level) return strings[6] ; 
  if ( G_LOG_LEVEL_DEBUG & level) return strings[7] ; 

  return NULL ;
}

/** 
 * Logging function for WBFMM. 
 * 
 * @param log_domain default log domain (see glib documentation);
 * @param log_level logging level (see glib documentation);
 * @param message logging message;
 * @param data array containing logging data, set by ::wbfmm_logging_init. 
 */

void wbfmm_logging_func(const gchar *log_domain,
		       GLogLevelFlags log_level,
		       const gchar *message,
		       gpointer data[])

{
  FILE *f = (FILE *)data[WBFMM_LOGGING_DATA_FID] ;
  gchar *p = (gchar *)data[WBFMM_LOGGING_DATA_PREFIX] ;
  GLogLevelFlags level = *(GLogLevelFlags *)data[WBFMM_LOGGING_DATA_LEVEL] ;
  gint (*exit_func)(void) = data[WBFMM_LOGGING_DATA_EXIT_FUNC] ;
  GTimer *timer = data[WBFMM_LOGGING_DATA_TIMER] ;
  
  if ( log_level > level ) return ;

  if ( timer == NULL ) 
    fprintf(f, "%s%s-%s: %s\n", p, 
	    G_LOG_DOMAIN, wbfmm_logging_string(log_level),
	    message) ;
  else
    fprintf(f, "%s%s-%s: %s [%lg]\n", p, 
	    G_LOG_DOMAIN, wbfmm_logging_string(log_level),
	    message, g_timer_elapsed(timer, NULL)) ;
    
  if ( log_level <= G_LOG_LEVEL_ERROR ) {
    if ( exit_func != NULL ) exit_func() ;
  }

  return ;
}

/** 
 * Initialize WBFMM logging
 * 
 * @param f file stream for messages;
 * @param p string to prepend to messages;
 * @param log_level maximum logging level to handle (see gts_log);
 * @param exit_func function to call if exiting on an error;
 * @param timed if TRUE time in seconds from start of logging is written.
 * 
 * @return 0 on success
 */

gint wbfmm_logging_init(FILE *f, gchar *p, 
		      GLogLevelFlags log_level,
		      gpointer exit_func, gboolean timed)

{
  static gpointer data[WBFMM_LOGGING_DATA_WIDTH] ;
  static GLogLevelFlags level ;
  GTimer *timer ;

  if ( f != NULL ) 
    data[WBFMM_LOGGING_DATA_FID] = f ;
  else
    data[WBFMM_LOGGING_DATA_FID] = stderr ;    
  if ( p != NULL ) 
    data[WBFMM_LOGGING_DATA_PREFIX] = g_strdup(p) ;
  else
    data[WBFMM_LOGGING_DATA_PREFIX] = g_strdup("") ;

  level = log_level ;
  data[WBFMM_LOGGING_DATA_LEVEL] = &level ;    

  if ( timed ) {
    timer = g_timer_new() ;
    g_timer_start(timer) ;
    data[WBFMM_LOGGING_DATA_TIMER] = timer ;
  } else {
    data[WBFMM_LOGGING_DATA_TIMER] = NULL ;
  }

  g_log_set_handler(G_LOG_DOMAIN, 
		    G_LOG_FLAG_RECURSION |
		    G_LOG_FLAG_FATAL |   
		    G_LOG_LEVEL_ERROR |
		    G_LOG_LEVEL_CRITICAL |
		    G_LOG_LEVEL_WARNING |
		    G_LOG_LEVEL_MESSAGE |
		    G_LOG_LEVEL_INFO |
		    G_LOG_LEVEL_DEBUG,
		    (GLogFunc)wbfmm_logging_func, data);

  return 0 ;
}

/**
 * @}
 * 
 */
