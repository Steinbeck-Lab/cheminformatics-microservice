// This file sets up an Axios instance for making API requests to the Natural Products API.
import axios, { AxiosInstance, InternalAxiosRequestConfig, AxiosResponse, AxiosError } from "axios";
import { toast } from "sonner";

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

    const status = axiosError.response?.status;

    // Network errors (no response) or 5xx server errors -- show global toast
    if (!axiosError.response || (status && status >= 500)) {
      toast.error("Service temporarily unavailable", {
        id: "network-error",
        description: "Please try again in a moment.",
        duration: 5000,
      });
    }

    // Rate limiting -- deduplicated toast
    if (status === 429) {
      toast.error("Too many requests", {
        id: "rate-limit",
        description: "Please wait before trying again.",
        duration: 5000,
      });
    }

    // 4xx errors are NOT toasted -- they propagate to per-component error handling
    return Promise.reject(error);
  }
);

export default api;
