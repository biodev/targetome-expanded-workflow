#!/usr/bin/env Rscript
#
# validate_pubmed_ids.R
# 
# Description: Validates a set of candidate PubMed IDs (PMIDs) by checking
#              whether they correspond to actual records in the PubMed database.
#
# Author: Anthropic Claude Opus 4.5
# Dependencies: None (uses only base R)
#
# Usage:
#   - Source the script and call validate_pubmed_ids() with a vector of PMIDs



#' Parse Simple JSON Response (Base R Implementation)
#'
#' A minimal JSON parser for the specific NCBI esummary response format.
#' Handles the essential structure needed for PMID validation.
#'
#' @param json_text Character string containing JSON
#' @return A list with uids and error_pmids
#' @keywords internal
parse_esummary_json <- function(json_text) {
  # Extract UIDs array from the response
  uids_match <- regmatches(json_text, 
                           regexpr('"uids"\\s*:\\s*\\[([^\\]]*)]', json_text))
  
  uids <- character(0)
  if (length(uids_match) > 0 && nzchar(uids_match)) {
    # Extract the array contents
    uids_content <- sub('.*\\[([^\\]]*)\\].*', '\\1', uids_match)
    # Parse individual UIDs (they're quoted strings)
    uids <- gsub('"', '', unlist(strsplit(uids_content, ",")))
    uids <- trimws(uids)
    uids <- uids[nzchar(uids)]
  }
  
  # Check for errors for each PMID
  # Pattern: "12345":{"uid":"12345","error":"..."} indicates error
  error_pattern <- '"([0-9]+)"\\s*:\\s*\\{[^}]*"error"\\s*:'
  error_matches <- gregexpr(error_pattern, json_text, perl = TRUE)
  error_pmids <- character(0)
  
  if (error_matches[[1]][1] != -1) {
    matches <- regmatches(json_text, error_matches)[[1]]
    error_pmids <- sub('"([0-9]+)".*', '\\1', matches)
  }
  
  return(list(uids = uids, error_pmids = error_pmids))
}


#' Build URL with Query Parameters
#'
#' @param base_url The base URL
#' @param params Named list of query parameters
#' @return Complete URL string
#' @keywords internal
build_url <- function(base_url, params) {
  if (length(params) == 0) return(base_url)
  
  query_parts <- sapply(names(params), function(name) {
    paste0(URLencode(name, reserved = TRUE), "=", 
           URLencode(as.character(params[[name]]), reserved = TRUE))
  })
  
  paste0(base_url, "?", paste(query_parts, collapse = "&"))
}


#' Make HTTP GET Request (Base R Implementation)
#'
#' @param url The URL to request
#' @param timeout_seconds Request timeout
#' @return List with status and content, or NULL on error
#' @keywords internal
http_get <- function(url, timeout_seconds = 30) {
  tryCatch({
    # Use base R url() connection
    con <- url(url, open = "rb")
    on.exit(close(con))
    
    # Read the response
    response_raw <- readLines(con, warn = FALSE)
    response_text <- paste(response_raw, collapse = "\n")
    
    return(list(
      status = 200,  # If we get here, assume success
      content = response_text
    ))
    
  }, error = function(e) {
    # Check if it's an HTTP error
    if (grepl("404|400|500|503", e$message)) {
      return(list(status = as.integer(sub(".*([0-9]{3}).*", "\\1", e$message)),
                  content = NULL))
    }
    warning(sprintf("HTTP request failed: %s", e$message))
    return(NULL)
  })
}


#' Validate Multiple PubMed IDs
#'
#' Validates a vector of PubMed IDs by querying the NCBI E-utilities API.
#' Uses batch queries for efficiency.
#'
#' @param pmids A vector of PubMed IDs (character or numeric)
#' @param api_key Optional NCBI API key for higher rate limits (default: NULL)
#'                Register at: https://www.ncbi.nlm.nih.gov/account/settings/
#' @param batch_size Number of PMIDs to query per batch (default: 200, max: 500)
#' @param delay_seconds Delay between batch requests to respect rate limits (default: 0.34)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A data.frame with columns:
#'   - pmid: The input PubMed ID
#'   - is_valid: TRUE if valid, FALSE if invalid, NA if error occurred
#'
#' @examples
#' # Validate a few PMIDs
#' results <- validate_pubmed_ids(c("12345678", "99999999", "33577776"))
#' print(results)
#'
#' # With an API key for faster processing
#' results <- validate_pubmed_ids(pmids, api_key = "your_api_key_here")
#'
#' @export
validate_pubmed_ids <- function(pmids, 
                                 api_key = NULL, 
                                 batch_size = 200,
                                 delay_seconds = 0.34,
                                 verbose = TRUE) {
  
  # Input validation
  if (length(pmids) == 0) {
    warning("No PMIDs provided")
    return(data.frame(pmid = character(0), is_valid = logical(0),
                      stringsAsFactors = FALSE))
  }
  
  # Sanitize input
  pmids <- as.character(pmids)
  pmids <- trimws(pmids)
  original_pmids <- pmids
  unique_pmids <- unique(pmids)
  
  if (verbose && length(unique_pmids) < length(pmids)) {
    message(sprintf("Note: %d duplicate PMIDs removed, validating %d unique PMIDs",
                    length(pmids) - length(unique_pmids), length(unique_pmids)))
  }
  
  # Pre-filter obviously invalid PMIDs (non-numeric or non-positive)
  is_numeric <- grepl("^[0-9]+$", unique_pmids)
  is_positive <- sapply(unique_pmids, function(x) {
    if (grepl("^[0-9]+$", x)) {
      # Handle very large numbers that might overflow
      tryCatch({
        as.numeric(x) > 0
      }, warning = function(w) TRUE)  # Assume valid format if parsing issues
    } else {
      FALSE
    }
  })
  format_valid <- is_numeric & is_positive
  pmids_to_check <- unique_pmids[format_valid]
  
  # Initialize results
  results <- rep(NA, length(unique_pmids))
  names(results) <- unique_pmids
  results[!format_valid] <- FALSE  # Mark format-invalid as FALSE
  
  if (verbose && sum(!format_valid) > 0) {
    message(sprintf("Note: %d PMIDs have invalid format (non-numeric or non-positive)",
                    sum(!format_valid)))
  }
  
  if (length(pmids_to_check) == 0) {
    if (verbose) message("No valid-format PMIDs to check against PubMed")
    result_df <- data.frame(
      pmid = original_pmids,
      is_valid = results[original_pmids],
      stringsAsFactors = FALSE
    )
    return(result_df)
  }
  
  # Batch processing
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
  n_batches <- ceiling(length(pmids_to_check) / batch_size)
  
  if (verbose) {
    message(sprintf("Validating %d PMIDs in %d batch(es)...", 
                    length(pmids_to_check), n_batches))
  }
  
  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(pmids_to_check))
    batch_pmids <- pmids_to_check[start_idx:end_idx]
    
    if (verbose) {
      message(sprintf("  Processing batch %d/%d (%d PMIDs)...", 
                      i, n_batches, length(batch_pmids)))
    }
    
    # Build query parameters
    query_params <- list(
      db = "pubmed",
      id = paste(batch_pmids, collapse = ","),
      retmode = "json"
    )
    
    if (!is.null(api_key) && nzchar(api_key)) {
      query_params$api_key <- api_key
    }
    
    # Make the request
    url <- build_url(base_url, query_params)
    
    tryCatch({
      response <- http_get(url)
      
      if (is.null(response)) {
        warning(sprintf("Request failed for batch %d", i))
        next
      }
      
      if (response$status != 200) {
        warning(sprintf("HTTP error %d for batch %d", response$status, i))
        next
      }
      
      # Parse the response
      parsed <- parse_esummary_json(response$content)
      
      # Mark valid PMIDs
      for (pmid in batch_pmids) {
        if (pmid %in% parsed$uids) {
          if (pmid %in% parsed$error_pmids) {
            results[pmid] <- FALSE
          } else {
            results[pmid] <- TRUE
          }
        } else {
          # Not in the returned UIDs - could be invalid or error
          # Check if there's an explicit error for this PMID
          if (pmid %in% parsed$error_pmids) {
            results[pmid] <- FALSE
          } else {
            # Double check: might not be in UIDs list but still have data
            # Search for the PMID's data in response
            pmid_pattern <- sprintf('"%s"\\s*:\\s*\\{[^}]*"uid"\\s*:\\s*"%s"', pmid, pmid)
            if (grepl(pmid_pattern, response$content)) {
              results[pmid] <- TRUE
            } else {
              results[pmid] <- FALSE
            }
          }
        }
      }
      
    }, error = function(e) {
      warning(sprintf("Error processing batch %d: %s", i, e$message))
    })
    
    # Rate limiting delay between batches
    if (i < n_batches) {
      Sys.sleep(delay_seconds)
    }
  }
  
  
  
  # Construct final results data frame maintaining original order and duplicates
  result_df <- data.frame(
    pmid = original_pmids,
    is_valid = results[original_pmids],
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- NULL
  
  if (verbose) {
    n_valid <- sum(result_df$is_valid == TRUE, na.rm = TRUE)
    n_invalid <- sum(result_df$is_valid == FALSE, na.rm = TRUE)
    n_error <- sum(is.na(result_df$is_valid))
    message(sprintf("Validation complete: %d valid, %d invalid, %d errors",
                    n_valid, n_invalid, n_error))
  }
  
  return(result_df)
}


#' Validate a Single PubMed ID
#'
#' Convenience function to check a single PMID.
#'
#' @param pmid A single PubMed ID (character or numeric)
#' @param api_key Optional NCBI API key
#' @return Logical: TRUE if valid, FALSE if invalid, NA on error
#' @export
is_valid_pmid <- function(pmid, api_key = NULL) {
  result <- validate_pubmed_ids(pmid, api_key = api_key, verbose = FALSE)
  return(result$is_valid[1])
}


#' Print Summary of Validation Results
#'
#' @param results Data frame returned by validate_pubmed_ids()
#' @export
summarize_validation <- function(results) {
  cat("\n=== PubMed ID Validation Summary ===\n")
  cat(sprintf("Total PMIDs checked:    %d\n", nrow(results)))
  cat(sprintf("Valid PMIDs:            %d\n", sum(results$is_valid == TRUE, na.rm = TRUE)))
  cat(sprintf("Invalid PMIDs:          %d\n", sum(results$is_valid == FALSE, na.rm = TRUE)))
  cat(sprintf("Errors (check failed):  %d\n", sum(is.na(results$is_valid))))
  
  invalid_pmids <- results$pmid[results$is_valid == FALSE]
  if (length(invalid_pmids) > 0 && length(invalid_pmids) <= 20) {
    cat(sprintf("\nInvalid PMIDs: %s\n", paste(invalid_pmids, collapse = ", ")))
  } else if (length(invalid_pmids) > 20) {
    cat(sprintf("\nFirst 20 invalid PMIDs: %s, ...\n", 
                paste(head(invalid_pmids, 20), collapse = ", ")))
  }
  
  valid_pmids <- results$pmid[results$is_valid == TRUE]
  if (length(valid_pmids) > 0 && length(valid_pmids) <= 20) {
    cat(sprintf("\nValid PMIDs: %s\n", paste(valid_pmids, collapse = ", ")))
  } else if (length(valid_pmids) > 20) {
    cat(sprintf("\nFirst 20 valid PMIDs: %s, ...\n", 
                paste(head(valid_pmids, 20), collapse = ", ")))
  }
  cat("\n")
  
  invisible(results)
}