// This file sets up an Axios instance for making API requests to the Natural Products API.
import axios from 'axios';

// Get the API URL from environment variables or use a default
const API_URL = process.env.REACT_APP_API_URL || 'https://api.naturalproducts.net/latest';

// Create an Axios instance with default config
const api = axios.create({
  baseURL: API_URL,
  headers: {
    'Content-Type': 'application/json',
  },
  timeout: 30000, // 30 seconds timeout
});

// Request interceptor for adding auth tokens or other headers
api.interceptors.request.use(
  (config) => {
    // Add any common request handling here
    return config;
  },
  (error) => {
    return Promise.reject(error);
  }
);

// Response interceptor for common error handling
api.interceptors.response.use(
  (response) => {
    return response;
  },
  (error) => {
    // Handle common error scenarios
    if (error.response) {
      // Server responded with a status other than 2xx
      console.error('API Error Response:', error.response.data);

      if (error.response.status === 429) {
        // Rate limiting
        console.warn('API rate limit reached');
      }
    } else if (error.request) {
      // Request was made but no response was received
      console.error('API No Response:', error.request);
    } else {
      // Something else happened while setting up the request
      console.error('API Request Error:', error.message);
    }

    return Promise.reject(error);
  }
);

export default api;