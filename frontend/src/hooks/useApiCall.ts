// Description: Custom hook for handling API calls with loading and error states
import { useState, useCallback } from "react";

/**
 * Custom hook for handling API calls with loading and error states
 * @param {Function} apiFunction - The API function to call
 * @param {Object} options - Options for the hook
 * @param {boolean} options.immediate - Whether to call the API immediately
 * @param {Array} options.initialArgs - Initial arguments for immediate API call
 * @param {Function} options.onSuccess - Callback function on successful API call
 * @param {Function} options.onError - Callback function on API call error
 * @returns {Object} - API call state and functions
 */
const useApiCall = (apiFunction, options = {}) => {
  const { immediate = false, initialArgs = [], onSuccess = null, onError = null } = options;

  const [isLoading, setIsLoading] = useState(immediate);
  const [error, setError] = useState(null);
  const [data, setData] = useState(null);
  const [lastCallTime, setLastCallTime] = useState(null);

  /**
   * Execute the API call
   * @param {...any} args - Arguments to pass to the API function
   * @returns {Promise<any>} - Promise resolving to the API result
   */
  const call = useCallback(
    async (...args) => {
      setIsLoading(true);
      setError(null);
      setLastCallTime(new Date());

      try {
        const result = await apiFunction(...args);
        setData(result);

        if (onSuccess) {
          onSuccess(result);
        }

        return result;
      } catch (err) {
        const errorMessage = err.message || "An error occurred during the API call";
        setError(errorMessage);

        if (onError) {
          onError(err);
        }

        throw err;
      } finally {
        setIsLoading(false);
      }
    },
    [apiFunction, onSuccess, onError]
  );

  /**
   * Reset all state
   */
  const reset = useCallback(() => {
    setIsLoading(false);
    setError(null);
    setData(null);
    setLastCallTime(null);
  }, []);

  // Make immediate API call if requested
  useState(() => {
    if (immediate) {
      call(...initialArgs).catch(() => {
        // Error is already handled and stored in state
      });
    }
  });

  return {
    call,
    data,
    isLoading,
    error,
    lastCallTime,
    reset,
  };
};

/**
 * Create a reusable API call function with built-in loading and error handling
 * @param {Function} apiFunction - The API function to wrap
 * @returns {Function} - Wrapped function that returns loading, error, and data
 */
export const createApiCall = (apiFunction) => {
  return async (...args) => {
    const [setLoading, setError, setData] = args.slice(-3);
    const apiArgs = args.slice(0, -3);

    try {
      setLoading(true);
      setError(null);

      const result = await apiFunction(...apiArgs);
      setData(result);
      return result;
    } catch (err) {
      const errorMessage = err.message || "An error occurred during the API call";
      setError(errorMessage);
      throw err;
    } finally {
      setLoading(false);
    }
  };
};

export default useApiCall;
