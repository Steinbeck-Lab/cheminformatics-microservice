// This file sets up an Axios instance for making API requests to the Natural Products API.
import axios, { AxiosInstance, InternalAxiosRequestConfig, AxiosResponse, AxiosError } from "axios";

// Get the API URL from environment variables or use a default
const API_URL: string =
  import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest";

// Create an Axios instance with default config
const api: AxiosInstance = axios.create({
  baseURL: API_URL,
  headers: {
    "Content-Type": "application/json",
  },
  timeout: 30000, // 30 seconds timeout
});

// Request interceptor for adding auth tokens or other headers
api.interceptors.request.use(
  (config: InternalAxiosRequestConfig) => {
    // Add any common request handling here
    return config;
  },
  (error: unknown) => {
    return Promise.reject(error);
  }
);

// Response interceptor for common error handling
api.interceptors.response.use(
  (response: AxiosResponse) => {
    return response;
  },
  (error: unknown) => {
    const axiosError = error as AxiosError;
    if (import.meta.env.DEV) {
      if (axiosError.response) {
        console.error("API Error Response:", axiosError.response.status);
      } else if (axiosError.request) {
        console.error("API No Response");
      } else if (error instanceof Error) {
        console.error("API Request Error:", error.message);
      }
    }

    if (axiosError.response && axiosError.response.status === 429) {
      console.warn("API rate limit reached");
    }

    return Promise.reject(error);
  }
);

export default api;
