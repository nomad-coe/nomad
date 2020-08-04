module.exports = {
    "extends": [
        "standard",
        "plugin:react/recommended"
    ],
    "parser": "babel-eslint",
    "parserOptions": {
        "ecmaFeatures": {
            "jsx": true,
            "modules": true
        }
    },
    "globals": {
        "fetch": false,
        "browser": true
    },
    "plugins": [
        "react", "react-hooks"
    ],
    "rules": {
        "space-before-function-paren": ["error", "never"],
        "camelcase": [0],
        "react-hooks/rules-of-hooks": "error",
        "react-hooks/exhaustive-deps": "warn",
        "react/display-name": ["error", "never"]
    },
    "settings": {
        "react": {
            "version": "detect"
        }
    }
}