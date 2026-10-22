import React from "react";
import ReactDOM from "react-dom/client";
import { MotionConfig } from "motion/react";
// Display face for big headlines only (self-hosted); UI text stays Inter.
import "@fontsource-variable/antonio";
import "./styles/tailwind.css";
import App from "./App";

const root = ReactDOM.createRoot(document.getElementById("root")!);
root.render(
  <React.StrictMode>
    {/* Honour the OS "reduce motion" setting: motion drops transform and layout animation */}
    <MotionConfig reducedMotion="user">
      <App />
    </MotionConfig>
  </React.StrictMode>
);
