import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { Link } from "react-router-dom";

// Animation Variants
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.7, ease: "easeOut" } },
};

const headerVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.7, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};

const contentVariants = {
  hidden: { opacity: 0, y: 30 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.7, delay: 0.3, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};

const PrivacyPolicy = () => {
  const [currentYear] = useState(new Date().getFullYear());
  const [isLoaded, setIsLoaded] = useState(false);

  useEffect(() => {
    // Add a slight delay to ensure smooth animations after component mounts
    const timer = setTimeout(() => {
      setIsLoaded(true);
    }, 100);
    return () => clearTimeout(timer);
  }, []);

  return (
    <motion.div
      className="flex flex-col min-h-screen bg-slate-50 dark:bg-gray-900 text-slate-800 dark:text-slate-100 font-sans overflow-x-hidden"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* Background Effects */}
      <div className="absolute inset-0 -z-20 overflow-hidden pointer-events-none">
        <div className="absolute inset-0 dark:bg-gradient-to-br dark:from-slate-900 dark:via-indigo-950/50 dark:to-slate-950 opacity-0 dark:opacity-100 transition-opacity duration-500"></div>
        <div className="absolute inset-0 bg-gradient-to-br from-blue-50/80 via-white to-indigo-50/80 opacity-100 dark:opacity-0 transition-opacity duration-500"></div>
        <div className="absolute inset-0 w-full h-full bg-[radial-gradient(#e5e7eb_1px,transparent_1px)] dark:bg-[radial-gradient(#1e293b30_1px,transparent_1px)] [background-size:20px_20px] opacity-50 dark:opacity-30"></div>
      </div>

      {/* Content Area */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 z-10 flex-grow">
        <div className="w-full max-w-4xl mx-auto">
          {/* Page Header */}
          <motion.div
            className="mb-12 text-center"
            variants={headerVariants}
            initial="hidden"
            animate={isLoaded ? "visible" : "hidden"}
          >
            <h1 className="text-4xl md:text-5xl font-bold mb-6">
              <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-600 to-indigo-600 dark:from-blue-400 dark:to-indigo-400">
                Privacy Policy
              </span>
            </h1>
            <p className="text-lg text-slate-600 dark:text-slate-300">
              Effective date: April 7, {currentYear}
            </p>
          </motion.div>

          {/* Main Content */}
          <motion.div
            className="bg-white dark:bg-slate-800 rounded-xl shadow-lg border border-slate-200 dark:border-slate-700 overflow-hidden"
            variants={contentVariants}
            initial="hidden"
            animate={isLoaded ? "visible" : "hidden"}
          >
            <div className="p-8 md:p-10">
              <div className="prose prose-slate dark:prose-invert max-w-none">
                <p className="text-slate-700 dark:text-slate-300">
                  Cheminformatics Microservice ("us", "we", or "our") operates the
                  https://api.naturalproducts.net/ website and related applications (hereinafter
                  referred to as the "Service").
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  This page informs you of our policies regarding the collection, use, and
                  disclosure of personal data when you use our Service and the choices you have
                  associated with that data.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  We use your data to provide and improve the Service. By using the Service, you
                  agree to the collection and use of information in accordance with this policy.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">Definitions</h2>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>
                    <strong>Service:</strong> Service refers to the https://api.naturalproducts.net/
                    website and related applications operated by the Steinbeck group at Friedrich
                    Schiller University Jena.
                  </li>
                  <li>
                    <strong>Personal Data:</strong> Personal Data means data about a living
                    individual who can be identified from those data (or from those and other
                    information either in our possession or likely to come into our possession).
                  </li>
                  <li>
                    <strong>Usage Data:</strong> Usage Data is data collected automatically either
                    generated by the use of the Service or from the Service infrastructure itself
                    (for example, the duration of a page visit).
                  </li>
                  <li>
                    <strong>Cookies:</strong> Cookies are small files stored on your device
                    (computer or mobile device).
                  </li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">
                  Information Collection and Use
                </h2>

                <p className="text-slate-700 dark:text-slate-300">
                  We collect limited information for the purposes of monitoring service usage and
                  improving our Service. Our approach prioritizes privacy, and we collect only what
                  is necessary for these purposes.
                </p>

                <h3 className="text-slate-700 dark:text-slate-300">Types of Data Collected</h3>

                <h4 className="text-slate-700 dark:text-slate-300">Usage Data</h4>

                <p className="text-slate-700 dark:text-slate-300">
                  We collect anonymized usage data to understand how our Service is accessed and
                  used. This information helps us improve our service and ensure it meets user
                  needs. This Usage Data may include:
                </p>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>IP address (anonymized)</li>
                  <li>Browser type and version</li>
                  <li>Pages of our Service that you visit</li>
                  <li>Time and date of your visit</li>
                  <li>Time spent on pages</li>
                </ul>

                <p className="text-slate-700 dark:text-slate-300">
                  For searches conducted on our platform, we do not track individual search history
                  on our servers. When using our chemical structure search features, search history
                  may be stored locally on your device for convenience, but this data is not
                  transmitted to our servers. Additionally, our system does not retain information
                  related to search queries that return zero results.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  We may collect and analyze aggregated search keywords to improve performance and
                  optimize the user experience. These anonymized keywords, along with related logs
                  and reports, are automatically deleted after 365 days.
                </p>

                <h4 className="text-slate-700 dark:text-slate-300">Tracking & Cookies Data</h4>

                <p className="text-slate-700 dark:text-slate-300">
                  We use cookies and similar tracking technologies to track activity on our Service
                  and to improve user experience. Cookies are files with a small amount of data
                  which may include an anonymous unique identifier. You can instruct your browser to
                  refuse all cookies or to indicate when a cookie is being sent. However, if you do
                  not accept cookies, you may not be able to use some portions of our Service.
                </p>

                <p className="text-slate-700 dark:text-slate-300">Examples of Cookies we use:</p>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>
                    <strong>Session Cookies:</strong> We use Session Cookies to operate our Service.
                  </li>
                  <li>
                    <strong>Preference Cookies:</strong> We use Preference Cookies to remember your
                    preferences and various settings.
                  </li>
                  <li>
                    <strong>Security Cookies:</strong> We use Security Cookies for security
                    purposes.
                  </li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">Analytics</h2>

                <p className="text-slate-700 dark:text-slate-300">
                  We use the following analytics service to monitor and analyze the use of our
                  Service:
                </p>

                <h3 className="text-slate-700 dark:text-slate-300">Matomo Analytics</h3>

                <p className="text-slate-700 dark:text-slate-300">
                  Matomo is used to analyze visitor behavior to identify potential improvements and
                  understand how users interact with our Service. All data is processed with privacy
                  in mind, and we employ the following privacy-protective measures:
                </p>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>IP addresses are anonymized</li>
                  <li>Data is stored on our own servers, not shared with third parties</li>
                  <li>We honor Do Not Track preferences</li>
                  <li>Analytics data is automatically deleted after 365 days</li>
                </ul>

                <p className="text-slate-700 dark:text-slate-300">
                  The processing of this data is based on our legitimate interest to improve our
                  Service and user experience.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">Use of Data</h2>

                <p className="text-slate-700 dark:text-slate-300">
                  We use the collected data for various purposes:
                </p>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>To provide and maintain our Service</li>
                  <li>To notify you about changes to our Service</li>
                  <li>
                    To allow you to participate in interactive features of our Service when you
                    choose to do so
                  </li>
                  <li>To provide support</li>
                  <li>
                    To gather analysis or valuable information so that we can improve our Service
                  </li>
                  <li>To monitor the usage of our Service</li>
                  <li>To detect, prevent, and address technical issues</li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">
                  Legal Basis for Processing Personal Data under the General Data Protection
                  Regulation (GDPR)
                </h2>

                <p className="text-slate-700 dark:text-slate-300">
                  If you are from the European Economic Area (EEA), our legal basis for collecting
                  and using the personal information described in this Privacy Policy depends on the
                  data we collect and the specific context in which we collect it.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  We may process your data because:
                </p>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>We need to perform a contract with you</li>
                  <li>You have given us permission to do so</li>
                  <li>
                    The processing is in our legitimate interests and is not overridden by your
                    rights
                  </li>
                  <li>To comply with the law</li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">Retention of Data</h2>

                <p className="text-slate-700 dark:text-slate-300">
                  We retain personal data only for as long as is necessary for the purposes set out
                  in this Privacy Policy. We retain and use your data to the extent necessary to
                  comply with our legal obligations, resolve disputes, and enforce our policies.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  Usage data is generally retained for a shorter period. We automatically delete
                  analytics data after 365 days.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">
                  Your Data Protection Rights under the General Data Protection Regulation (GDPR)
                </h2>

                <p className="text-slate-700 dark:text-slate-300">
                  If you are a resident of the European Economic Area (EEA), you have certain data
                  protection rights. We aim to take reasonable steps to allow you to correct, amend,
                  delete, or limit the use of your Personal Data.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  You have the following data protection rights:
                </p>

                <ul className="text-slate-700 dark:text-slate-300">
                  <li>The right to access, update or delete the information we have on you</li>
                  <li>The right of rectification</li>
                  <li>The right to object</li>
                  <li>The right of restriction</li>
                  <li>The right to data portability</li>
                  <li>The right to withdraw consent</li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">Links to Other Sites</h2>

                <p className="text-slate-700 dark:text-slate-300">
                  Our Service may contain links to other sites that are not operated by us. If you
                  click a third-party link, you will be directed to that third party's site. We
                  strongly advise you to review the Privacy Policy of every site you visit.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  We have no control over and assume no responsibility for the content, privacy
                  policies, or practices of any third-party sites or services.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">Children's Privacy</h2>

                <p className="text-slate-700 dark:text-slate-300">
                  Our Service does not address anyone under the age of 18 ("Children"). We do not
                  knowingly collect personally identifiable information from anyone under the age of
                  18.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">
                  Changes to This Privacy Policy
                </h2>

                <p className="text-slate-700 dark:text-slate-300">
                  We may update our Privacy Policy from time to time. We will notify you of any
                  changes by posting the new Privacy Policy on this page and updating the "effective
                  date" at the top of this Privacy Policy.
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  You are advised to review this Privacy Policy periodically for any changes.
                  Changes to this Privacy Policy are effective when they are posted on this page.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">Contact Us</h2>

                <p className="text-slate-700 dark:text-slate-300">
                  If you have any questions about this Privacy Policy, please contact the Steinbeck
                  group at the Friedrich Schiller University Jena:
                </p>

                <p className="text-slate-700 dark:text-slate-300">
                  <a
                    href="https://cheminf.uni-jena.de"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline"
                  >
                    https://cheminf.uni-jena.de
                  </a>
                </p>

                <div className="mt-8 pt-6 border-t border-slate-200 dark:border-slate-700">
                  <Link to="/terms" className="text-blue-600 dark:text-blue-400 hover:underline">
                    View our Terms of Service
                  </Link>
                </div>
              </div>
            </div>
          </motion.div>
        </div>
      </div>
    </motion.div>
  );
};

export default PrivacyPolicy;
