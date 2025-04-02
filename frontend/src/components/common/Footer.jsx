// Description: Footer component with enhanced animations, improved layout, and additional features for better user experience.
import React, { useEffect, useRef } from 'react';
import { FaGithub, FaBook, FaFlask, FaUniversity, FaCoffee, FaCode } from 'react-icons/fa';
import { motion, useAnimation, useInView } from 'framer-motion';
import { useAppContext } from '../../context/AppContext';
import { Link } from 'react-router-dom';

// Enhanced animation variants
const footerVariant = {
  hidden: { opacity: 0, y: 30 },
  visible: {
    opacity: 1,
    y: 0,
    transition: {
      duration: 0.8,
      ease: [0.25, 0.1, 0.25, 1],
      when: "beforeChildren"
    }
  }
};

// Improved stagger for container elements
const staggerContainer = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: {
      staggerChildren: 0.1,
      delayChildren: 0.2,
      ease: "easeOut"
    }
  }
};

// Enhanced individual item animation
const itemVariant = {
  hidden: { opacity: 0, y: 15, scale: 0.95 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: {
      duration: 0.6,
      ease: [0.17, 0.67, 0.83, 0.67]
    }
  }
};


// Simplified Particle effect component
const Particles = () => {
  const particlesRef = useRef(null);

  useEffect(() => {
    const canvas = particlesRef.current;
    const ctx = canvas.getContext('2d');
    let animationFrameId;

    canvas.width = canvas.offsetWidth;
    canvas.height = canvas.offsetHeight;

    const particles = [];
    const particleCount = 12; // Reduced count

    for (let i = 0; i < particleCount; i++) {
      particles.push({
        x: Math.random() * canvas.width,
        y: Math.random() * canvas.height,
        radius: Math.random() * 2 + 0.5, // Smaller particles
        color: `rgba(59, 130, 246, ${Math.random() * 0.4 + 0.1})`,
        speedX: Math.random() * 0.6 - 0.3, // Slower movement
        speedY: Math.random() * 0.6 - 0.3
      });
    }

    const animate = () => {
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      particles.forEach(particle => {
        particle.x += particle.speedX;
        particle.y += particle.speedY;

        if (particle.x < 0 || particle.x > canvas.width) particle.speedX *= -1;
        if (particle.y < 0 || particle.y > canvas.height) particle.speedY *= -1;

        ctx.beginPath();
        ctx.arc(particle.x, particle.y, particle.radius, 0, Math.PI * 2);
        ctx.fillStyle = particle.color;
        ctx.fill();
      });

      animationFrameId = requestAnimationFrame(animate);
    };

    animate();

    const handleResize = () => {
      canvas.width = canvas.offsetWidth;
      canvas.height = canvas.offsetHeight;
    };

    window.addEventListener('resize', handleResize);

    return () => {
      window.removeEventListener('resize', handleResize);
      cancelAnimationFrame(animationFrameId);
    };
  }, []);

  return (
    <canvas
      ref={particlesRef}
      className="absolute inset-0 w-full h-full opacity-30 pointer-events-none z-0"
    />
  );
};

// Main Footer Component
const Footer = () => {
  const { isDarkMode } = useAppContext();
  const currentYear = new Date().getFullYear();
  const footerRef = useRef(null);
  const isInView = useInView(footerRef, { once: false, amount: 0.1 });
  const controls = useAnimation();

  useEffect(() => {
    if (isInView) {
      controls.start("visible");
    }
  }, [controls, isInView]);

  // Enhanced logo URLs with animated versions
  const logoDark = "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";
  const logoLight = "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";

  // Enhanced Resource Links with additional icons & descriptions - more compact
  const resources = [
    {
      href: "https://api.naturalproducts.net/latest/docs",
      icon: FaBook,
      title: "API Docs",
      desc: "Explore endpoints",
      color: "from-blue-500 to-indigo-500"
    },
    {
      href: "https://github.com/Steinbeck-Lab/cheminformatics-microservice",
      icon: FaGithub,
      title: "GitHub",
      desc: "View code",
      color: "from-purple-500 to-pink-500"
    },
    {
      href: "https://docs.api.naturalproducts.net",
      icon: FaCode,
      title: "Guides",
      desc: "Integration help",
      color: "from-green-500 to-teal-500"
    },
    {
      href: "https://cheminf.uni-jena.de",
      icon: FaUniversity,
      title: "Steinbeck Lab",
      desc: "Research group",
      color: "from-amber-500 to-orange-500"
    }
  ];

  return (
    <motion.footer
      ref={footerRef}
      className="relative overflow-hidden bg-gradient-to-b from-slate-100 via-white to-slate-100 dark:from-gray-950 dark:via-slate-900 dark:to-gray-950 py-6 mt-auto border-t border-slate-200 dark:border-slate-800" // Reduced padding
      variants={footerVariant}
      initial="hidden"
      animate={controls}
    >
      {/* Animated background particles */}
      <Particles />

      {/* Glass morphism overlay for depth */}
      <div className="absolute inset-0 bg-white/30 dark:bg-black/20 backdrop-blur-[2px] z-0"></div>

      {/* Decorative molecule SVG in background - smaller and more transparent */}
      <div className="absolute -right-16 -bottom-16 opacity-3 dark:opacity-8 pointer-events-none scale-75">
        <svg width="200" height="200" viewBox="0 0 200 200" fill="none" xmlns="http://www.w3.org/2000/svg">
          <circle cx="100" cy="100" r="10" fill="currentColor" className="text-blue-600" />
          <circle cx="150" cy="70" r="8" fill="currentColor" className="text-purple-600" />
          <circle cx="60" cy="130" r="8" fill="currentColor" className="text-green-600" />
          <circle cx="130" cy="140" r="8" fill="currentColor" className="text-amber-600" />
          <line x1="100" y1="100" x2="150" y2="70" stroke="currentColor" strokeWidth="2" className="text-gray-400" />
          <line x1="100" y1="100" x2="60" y2="130" stroke="currentColor" strokeWidth="2" className="text-gray-400" />
          <line x1="100" y1="100" x2="130" y2="140" stroke="currentColor" strokeWidth="2" className="text-gray-400" />
        </svg>
      </div>

      {/* Main content - more compact */}
      <motion.div
        className="relative max-w-7xl mx-auto px-4 lg:px-6 z-10"
        variants={staggerContainer}
      >
        {/* Top Section with Logo and Resources - more compact layout */}
        <div className="flex flex-col md:flex-row justify-between items-center md:items-start gap-6">
          {/* Logo and Description - more compact */}
          <motion.div
            className="md:w-1/3 mb-2 md:mb-0 text-center md:text-left"
            variants={itemVariant}
          >
            {/* Animated Logo */}
            <Link to="/" className="inline-flex items-center group" aria-label="Homepage">
              <motion.div
                className="relative overflow-hidden rounded-lg"
                whileHover={{ scale: 1.05 }}
                transition={{ type: 'spring', stiffness: 300, damping: 10 }}
              >
                {/* Glowing animated border */}
                <motion.div
                  className="absolute inset-0 bg-gradient-to-r from-blue-500 via-purple-500 to-blue-500 opacity-75 rounded-lg"
                  animate={{
                    backgroundPosition: ['0% 50%', '100% 50%', '0% 50%'],
                  }}
                  transition={{
                    duration: 5,
                    repeat: Infinity,
                    ease: "linear"
                  }}
                  style={{
                    backgroundSize: '200% 200%',
                  }}
                />

                <motion.img
                  src={isDarkMode ? logoDark : logoLight}
                  alt="Cheminformatics Logo"
                  className="h-8 w-auto relative z-10 bg-white dark:bg-gray-900 p-1 rounded-md" // Smaller image
                  whileHover={{
                    rotate: [-3, 3, -3],
                    transition: { duration: 1.5, repeat: Infinity }
                  }}
                />
              </motion.div>

              <div className="ml-2">
                {/* Animated gradient text */}
                <motion.span
                  className="font-bold text-lg block leading-tight bg-clip-text text-transparent" // Smaller text
                  style={{
                    backgroundImage: 'linear-gradient(90deg, #3b82f6, #8b5cf6, #3b82f6)',
                    backgroundSize: '200% auto',
                  }}
                  animate={{
                    backgroundPosition: ['0% center', '200% center', '0% center'],
                  }}
                  transition={{
                    duration: 5,
                    repeat: Infinity,
                    ease: "linear",
                  }}
                >
                  Cheminformatics
                </motion.span>
                <motion.span
                  className="text-xs text-slate-500 dark:text-blue-300 group-hover:text-blue-600 dark:group-hover:text-blue-200 transition-colors leading-tight"
                  whileHover={{ scale: 1.05, x: 2 }}
                >
                  Microservices
                </motion.span>
              </div>
            </Link>

            {/* Brief description - more compact */}
            <p className="text-slate-700 dark:text-slate-300 text-xs leading-tight mt-1 ml-1 max-w-xs"> {/* Smaller text, tighter leading */}
              Modern interface for chemical data analysis, conversion, and visualization
              powered by the Cheminformatics Microservice API.
            </p>
          </motion.div>

          {/* Resources with Enhanced Icons - more compact */}
          <motion.div
            className="md:w-2/3 grid grid-cols-2 sm:grid-cols-4 gap-3"
            variants={itemVariant}
          >
            {resources.map((link, index) => (
              <motion.a
                key={index}
                href={link.href}
                target="_blank"
                rel="noopener noreferrer"
                // MODIFIED: Reduced padding from p-3 to p-1
                className="group flex flex-col justify-center items-start p-1 rounded-lg border border-slate-200 dark:border-slate-700 bg-white/50 dark:bg-slate-800/50 backdrop-blur-sm transition-all duration-300 hover:shadow-md hover:border-blue-300 dark:hover:border-blue-500"
                variants={itemVariant}
                whileHover={{
                  y: -3,
                  boxShadow: "0 10px 25px -5px rgba(59, 130, 246, 0.15)"
                }}
                whileTap={{ scale: 0.98 }}
              >
                {/* Gradient accent line */}
                {/* MODIFIED: Reduced margin-bottom from mb-2 to mb-1 */}
                <div className={`h-1 w-8 mb-1 rounded-full bg-gradient-to-r ${link.color}`}></div>

                {/* Icon and title in horizontal layout */}
                {/* MODIFIED: Reduced margin-bottom from mb-1 to mb-0.5 */}
                <div className="flex items-center gap-2 mb-0.5">
                  <motion.div
                    className="text-slate-700 dark:text-slate-300 group-hover:text-blue-600 dark:group-hover:text-blue-400 transition-colors"
                    whileHover={{ scale: 1.1, rotate: 5 }}
                  >
                    {/* Icon size kept the same, could be reduced if needed */}
                    <link.icon className="h-4 w-4" />
                  </motion.div>
                  {/* Font size kept the same, could be reduced if needed */}
                  <h4 className="font-medium text-slate-900 dark:text-white text-sm">{link.title}</h4>
                </div>

                {/* Description */}
                {/* Font size kept the same, could be reduced if needed */}
                <p className="text-xs text-slate-600 dark:text-slate-400">
                  {link.desc}
                </p>

                {/* Animated arrow that appears on hover */}
                {/* MODIFIED: Reduced margin-top from mt-1 to mt-0.5 and height from h-4 to h-3 */}
                <motion.div
                  className="mt-0.5 overflow-hidden h-3 w-full" // Reduced height here
                  initial={{ opacity: 0 }}
                  whileHover={{ opacity: 1 }}
                >
                  <motion.span
                    className="text-blue-600 dark:text-blue-400 text-xs flex items-center"
                    initial={{ x: -20, opacity: 0 }}
                    whileHover={{ x: 0, opacity: 1 }}
                    transition={{ duration: 0.3 }}
                  >
                    Explore →
                  </motion.span>
                </motion.div>
              </motion.a>
            ))}
          </motion.div>
        </div>

        {/* Enhanced Divider with Animated Icon - more compact */}
        <motion.div
          className="relative py-2 flex items-center justify-center mt-3" // Reduced padding
          variants={itemVariant}
        >
          <div className="flex-grow h-px bg-gradient-to-r from-transparent via-slate-300 dark:via-slate-600 to-transparent"></div>

          <div className="mx-4"> {/* Reduced margin */}
            <motion.div
              className="relative p-1 rounded-full" // Smaller padding
              animate={{
                boxShadow: [
                  '0 0 0 rgba(59, 130, 246, 0)',
                  '0 0 10px rgba(59, 130, 246, 0.5)',
                  '0 0 0 rgba(59, 130, 246, 0)'
                ]
              }}
              transition={{
                duration: 3,
                repeat: Infinity,
                ease: "easeInOut"
              }}
            >
              {/* Background glow */}
              <motion.div
                className="absolute inset-0 rounded-full bg-blue-500/20 dark:bg-blue-600/20 backdrop-blur-sm"
                animate={{
                  scale: [1, 1.2, 1],
                  opacity: [0.6, 0.8, 0.6]
                }}
                transition={{
                  duration: 3,
                  repeat: Infinity,
                  ease: "easeInOut"
                }}
              />

              {/* Rotating and pulsing flask icon - smaller */}
              <motion.div
                className="relative z-10 bg-white dark:bg-gray-800 p-1 rounded-full" // Smaller padding
                animate={{
                  rotate: [0, 10, -5, 0],
                  scale: [1, 1.1, 0.95, 1]
                }}
                transition={{
                  duration: 5,
                  repeat: Infinity,
                  ease: "easeInOut"
                }}
              >
                <FaFlask className="h-4 w-4 text-blue-600 dark:text-blue-400" /> {/* Smaller icon */}
              </motion.div>
            </motion.div>
          </div>

          <div className="flex-grow h-px bg-gradient-to-r from-transparent via-slate-300 dark:via-slate-600 to-transparent"></div>
        </motion.div>

        {/* Enhanced Copyright section - more compact */}
        <motion.div
          className="mt-2 text-center relative" // Reduced margin
          variants={itemVariant}
        >
          {/* Decorative gradient line - smaller */}
          <motion.div
            className="h-0.5 w-16 mx-auto mb-2 rounded-full bg-gradient-to-r from-blue-500 via-purple-500 to-blue-500" // Thinner line, reduced margin
            animate={{
              backgroundPosition: ['0% 50%', '100% 50%', '0% 50%'],
            }}
            transition={{
              duration: 5,
              repeat: Infinity,
              ease: "linear"
            }}
            style={{
              backgroundSize: '200% 200%',
            }}
          />

          <p className="text-slate-600 dark:text-slate-300 text-[10px] flex items-center justify-center gap-1"> {/* Smaller text, reduced gap */}
            &copy; {currentYear} Cheminformatics Microservice UI. Built with
            <motion.span
              className="inline-flex items-center"
              animate={{
                y: [-1, 1, -1],
                rotate: [0, 5, 0, -5, 0],
              }}
              transition={{
                y: {
                  repeat: Infinity,
                  duration: 1.5,
                  ease: "easeInOut"
                },
                rotate: {
                  repeat: Infinity,
                  duration: 2,
                  ease: "easeInOut",
                  delay: 0.5
                }
              }}
            >
              <FaCoffee className="text-amber-600 dark:text-amber-400 h-3 w-3" /> {/* Smaller icon */}
            </motion.span>
            and the Cheminformatics Microservice API.
          </p>

          <motion.p
            className="mt-1 text-[8px] text-slate-500 dark:text-slate-400" // Smaller text, reduced margin
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ delay: 0.5, duration: 0.5 }}
          >
            Developed by the Steinbeck Lab at Friedrich Schiller University Jena.
            Funded by the Deutsche Forschungsgemeinschaft (DFG).
          </motion.p>
        </motion.div>
      </motion.div>
    </motion.footer>
  );
};

export default Footer;